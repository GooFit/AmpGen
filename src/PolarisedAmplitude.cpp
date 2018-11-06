#include "AmpGen/PolarisedAmplitude.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Array.h"

using namespace AmpGen;

std::vector<int> polarisations( const std::string& name ){
  auto props = *ParticlePropertiesList::get( name );
  std::vector<int> rt( props.twoSpin() + 1 );

  if( props.isFermion() ) return {1,-1};    
  if( props.twoSpin() == 0 ) return {0};
  if( props.twoSpin() == 2 ) 
    return (name == "gamma0") ? std::vector<int>({-1,1}) : std::vector<int>({-1,0,1}) ;
  else return {0};
}

std::vector< std::vector< int > > polarisation_outer_product( const std::vector< std::vector< int > >& A, const std::vector< int >& B ){
  std::vector< std::vector< int > > rt; 
  for( auto& iA : A ){
    for( auto& iB : B ){
      rt.push_back( iA );
      rt.rbegin()->push_back( iB );
    }
  }
  return rt;
}

PolarisedAmplitude::PolarisedAmplitude( const EventType& type, AmpGen::MinuitParameterSet& mps, const std::string& prefix ) : 
  m_nCalls(0),
  m_events(nullptr),
  m_weightParam(nullptr),
  m_mps(&mps),
  m_weight(1),
  m_eventType(type)
{
  AmplitudeRules proto( mps ); 

  auto proto_amplitudes = proto.getMatchingRules( type, prefix , true );

  auto production_polarisations = polarisations( type.mother() ); 
  std::vector<std::vector<int>> all_polarisation_states ;
  for( auto& pol : production_polarisations ) all_polarisation_states.push_back( {pol} );

  for( unsigned int i = 0 ; i < type.size(); ++i ){
    all_polarisation_states = polarisation_outer_product( all_polarisation_states, polarisations( type[i] ) );
  }
  auto set_polarisation_state = []( auto& matrix_element, auto& polState ){
    auto fs = matrix_element.first.getFinalStateParticles();
    matrix_element.first.setPolarisationState( polState[0] );
    for( unsigned int i = 0 ; i < fs.size(); ++i ){
      fs[i]->setPolarisationState( polState[i+1] );
    }
    matrix_element.first.makeUniqueString(); 
  };
  bool debug = NamedParameter<bool>("PolarisedAmplitude::Debug", false );
  for( auto& matrix_element : proto_amplitudes ){
    Tensor thisExpression( std::vector<size_t>({all_polarisation_states.size()}) );
    int i = 0 ;
    DebugSymbols syms;  
    for( auto& polState : all_polarisation_states ){
      set_polarisation_state( matrix_element, polState );
      DEBUG( "Calculating matrix element for: [" << vectorToString( polState, ", " ) << "] " << matrix_element.first.uniqueString()  );
      thisExpression[ i++ ] = make_cse( matrix_element.first.getExpression(&syms) ); 
    }
    CompiledExpression< std::vector<complex_t> , const real_t*, const real_t* > expression( 
        TensorExpression( thisExpression), 
        "p"+std::to_string( FNV1a_hash( matrix_element.first.uniqueString() ) ) , 
        type.getEventFormat(),debug ? syms : DebugSymbols() ,&mps ); 
    INFO( matrix_element.first.uniqueString() ); 
    m_matrixElements.emplace_back( std::make_shared<Particle>(matrix_element.first), matrix_element.second, expression );
  } 
  for( auto& polState : all_polarisation_states){
    for( auto& matrix_element : proto_amplitudes ){
      set_polarisation_state( matrix_element, polState );
    }
  }
  m_polStates = all_polarisation_states; 
  for( auto& thing : m_matrixElements ) CompilerWrapper().compile( thing.pdf, "" );
  if( mps.find("Px") == nullptr ){
    WARNING("Polarisation parameters not defined, defaulting to (0,0,0)");
  }
  m_productionParameters.emplace_back( mps.addOrGet("Px",2,0,0 ) );
  m_productionParameters.emplace_back( mps.addOrGet("Py",2,0,0 ) );
  m_productionParameters.emplace_back( mps.addOrGet("Pz",2,0,0 ) ); 
  for( size_t i = 0 ; i < 6 ; ++i ) m_norms[i].resize( m_matrixElements.size(), m_matrixElements.size() );
}

std::vector<TransitionMatrix<std::vector<complex_t>>> PolarisedAmplitude::matrixElements() const
{
  return m_matrixElements;  
}

void   PolarisedAmplitude::prepare()
{
  bool isReady = 0 ; 
  auto dim = m_eventType.dim();
  std::vector<size_t> changedPdfIndices; 
  auto tStartEval = std::chrono::high_resolution_clock::now();
  for( size_t i = 0; i < m_matrixElements.size(); ++i ){
    auto tBeginCal = std::chrono::high_resolution_clock::now();
    auto& t = m_matrixElements[i];
    if( m_nCalls == 0 ) isReady &= t.pdf.isReady(); 
    t.pdf.prepare();
    if( m_nCalls != 0 && !t.pdf.hasExternalsChanged() ) continue; 
    if( t.addressData == 999 ){
      DEBUG("Registering expression: " << t.decayTree->uniqueString() << " = " << t.addressData );
      t.addressData = m_events->registerExpression( t.pdf , dim.first * dim.second );
    }
    m_events->updateCache( t.pdf, t.addressData );
    double timeNorm   = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tBeginCal ).count();
    t.pdf.resetExternals();
    changedPdfIndices.push_back(i);
    DEBUG("Updated cache for PDF: " << std::setw(55) << t.decayTree->uniqueString() << std::setw(3) << t.addressData << " " << "  t = " << timeNorm << " ms" );
  }
  if( !m_probExpression.isLinked() ) build_probunnormalised();
  m_probExpression.prepare(); 
  m_weight = m_weightParam == nullptr ? 1 : m_weightParam->mean();
  auto tEndEval        = std::chrono::high_resolution_clock::now();
  
  if( m_integralDispatch.sim != nullptr )
  {
    if( changedPdfIndices.size() != 0 ) calculateNorms(changedPdfIndices);
    double px = m_productionParameters[0];
    double py = m_productionParameters[1];
    double pz = m_productionParameters[2];
    std::array< complex_t, 6 > acc;
    for( size_t i = 0 ; i < m_matrixElements.size(); ++i )
    {
      for( size_t j = 0 ; j < m_matrixElements.size(); ++j )
      {
        complex_t c = m_matrixElements[i].coupling() * std::conj( m_matrixElements[j].coupling() );
        for( size_t k = 0 ; k < 6 ; ++k ) acc[k] += c * m_norms[k].get(i,j);
      }
    }
    complex_t z(px,-py);
    complex_t total = (1+pz)*( acc[0] + acc[1] )
      + (1-pz)*( acc[2] + acc[3] ) 
      + std::conj(z)*( acc[4] +acc[5] ) 
      + z * ( std::conj( acc[4] + acc[5] ) );
    m_norm = std::real(total);
    if( m_nCalls == 0 ) debug_norm();
  }
  auto tEndIntegral = std::chrono::high_resolution_clock::now();
  double timeNorm   = std::chrono::duration<double, std::milli>( tEndIntegral - tEndEval ).count();
  double timeEval   = std::chrono::duration<double, std::milli>( tEndEval     - tStartEval ).count();
  if( changedPdfIndices.size() != 0  ) INFO("Time to evaluate = " << timeEval << " ms; norm = " << timeNorm << " ms;  [pdfs = " << changedPdfIndices.size() << "]" );
  m_nCalls++;
}

void PolarisedAmplitude::debug_norm()
{
  double norm = 0;
  for( auto& evt : * m_integralDispatch.sim ){
    norm += evt.weight() * prob_unnormalised(evt) / evt.genPdf(); 
    if( prob_unnormalised(evt) > 1e6 ){
      WARNING("Event has very high probability: " << prob_unnormalised(evt) );
      evt.print();
    }
  }
  INFO("Average = " << m_norm );
  auto evt = (*m_integralDispatch.sim)[0];
  evt.print();
  INFO("Check one event: " << prob_unnormalised(evt) << " " << getValNoCache(evt) );
  INFO("Cross-checking normalisation: " << m_norm << " vs " << norm / m_integralDispatch.sim->norm() << " sample_norm = " << m_integralDispatch.sim->norm() );
  int nNorms = 0 ; 
  int nZeros = 0; 
  for( size_t i = 0 ; i < m_matrixElements.size(); ++i ){
    for( size_t j = 0 ; j < m_matrixElements.size(); ++j ){
      std::string norm_string = "[";
      for( size_t k = 0 ; k < 6 ; ++k ) {
        double re = std::real( m_norms[k].get(i,j) );
        double im = std::imag( m_norms[k].get(i,j) );
        norm_string += "(" + std::to_string(re) +", " + std::to_string(im) + ") , ";
        if( re < 1e-8 && im < 1e-8 ) nZeros ++;
      }
      norm_string += "]";
      nNorms += 6; 
      INFO( m_matrixElements[i].decayTree->uniqueString() << " " << m_matrixElements[j].decayTree->uniqueString() << " " << norm_string );
    }
  }
  INFO( "nNorms = " << nNorms << " nZeros = " <<nZeros );
}

void   PolarisedAmplitude::setEvents( AmpGen::EventList& events )
{ 
  reset();
  m_events = &events;
}

void   PolarisedAmplitude::setMC( AmpGen::EventList& events )
{
  m_nCalls = 0;
  m_integralDispatch.sim = &events;
}
size_t PolarisedAmplitude::size() const 
{ 
  auto dim = m_eventType.dim() ; 
  return dim.first * dim.second * m_matrixElements.size(); 
}
void   PolarisedAmplitude::reset( const bool& flag ){ m_nCalls = 0 ; }

void PolarisedAmplitude::build_probunnormalised()
{
  Expression prob = probExpression( transitionMatrix() , { Parameter("Px"), Parameter("Py"), Parameter("Pz") } );
  m_probExpression = CompiledExpression< double, const real_t*, const complex_t*>( 
      prob, "prob_unnormalised", std::map<std::string, size_t>(), {}, m_mps );
   CompilerWrapper().compile( m_probExpression );
} 
  
Tensor PolarisedAmplitude::transitionMatrix()
{ 
  auto dim = m_eventType.dim();
  auto size = dim.first * dim.second ; 
  std::vector< Expression > expressions( size, 0);
  for( auto& me : m_matrixElements ){
    auto coupling   = me.coupling.to_expression() ;
    auto cacheIndex = m_events->getCacheIndex( me.pdf ); 
    for( size_t i = 0 ; i < size ; ++i ){
      expressions[i] = expressions[i] + coupling * Parameter( "x1["+std::to_string(cacheIndex+i)+"]",0,true); 
    }
  } 
  Tensor T_matrix( expressions, { dim.first , dim.second } );
  T_matrix.st();
  return T_matrix; 
}

complex_t PolarisedAmplitude::TE( const Event& event, const size_t& x, const size_t& y )
{
  auto dim = m_eventType.dim();
  size_t size = dim.first * dim.second; 
  std::vector<complex_t> acc( size ); 
  for( auto& me : m_matrixElements ){
    auto coupling   = me.coupling();
    auto cacheIndex = m_events->getCacheIndex( me.pdf ); 
    for( size_t i = 0 ; i < size ; ++i ){
      acc[i] = acc[i] + coupling * event.getCache( cacheIndex + i );
    }
  }
  complex_t T; 
  for( size_t k = 0; k < dim.second; ++k )
    T += acc[ x*dim.second +k ] * std::conj( acc[y*dim.second +k ] );
  return T;
}


double PolarisedAmplitude::prob_unnormalised( const AmpGen::Event& evt ) const
{
  return m_probExpression( evt.getCachePtr(0) );
}

double PolarisedAmplitude::norm() const 
{
  return m_norm;
}

struct pclock {
  std::chrono::_V2::system_clock::time_point t_start; 
  std::chrono::_V2::system_clock::time_point t_end; 
  pclock() : t_start(   std::chrono::high_resolution_clock::now()) {}
  void stop(){ t_end = std::chrono::high_resolution_clock::now() ; }
  operator double() const { return std::chrono::duration<double, std::milli>( t_end - t_start ).count() ; } ; 
};

void PolarisedAmplitude::calculateNorms(const std::vector<size_t>& changedPdfIndices ) 
{
  size_t size_of = size() / m_matrixElements.size();
  pclock t_eval; 
  for( auto& i : changedPdfIndices ){
    m_integralDispatch.prepareExpression( m_matrixElements[i].pdf , size_of );
  }
  t_eval.stop();
  pclock t_integral; 
  size_t total_integrals = 0; 
  std::vector<bool> integral_has_changed( 6 * m_matrixElements.size() * m_matrixElements.size() );
  for( auto& i : changedPdfIndices ){
    auto ai = m_integralDispatch.sim->getCacheIndex( m_matrixElements[i].pdf );
    for(size_t j = 0; j < m_matrixElements.size(); ++j ){
      auto aj = m_integralDispatch.sim->getCacheIndex( m_matrixElements[j].pdf );
      m_integralDispatch.addIntegralKeyed(ai+0, aj+0, [i,j,this](const complex_t& val){this->m_norms[0].set(i,j,val); if( i != j ) this->m_norms[0].set(j,i,std::conj(val));});
      m_integralDispatch.addIntegralKeyed(ai+1, aj+1, [i,j,this](const complex_t& val){this->m_norms[1].set(i,j,val); if( i != j ) this->m_norms[1].set(j,i,std::conj(val));});
      m_integralDispatch.addIntegralKeyed(ai+2, aj+2, [i,j,this](const complex_t& val){this->m_norms[2].set(i,j,val); if( i != j ) this->m_norms[2].set(j,i,std::conj(val));});
      m_integralDispatch.addIntegralKeyed(ai+3, aj+3, [i,j,this](const complex_t& val){this->m_norms[3].set(i,j,val); if( i != j ) this->m_norms[3].set(j,i,std::conj(val));});
      m_integralDispatch.addIntegralKeyed(ai+0, aj+2, [i,j,this](const complex_t& val){this->m_norms[4].set(i,j,val); }); 
      m_integralDispatch.addIntegralKeyed(ai+1, aj+3, [i,j,this](const complex_t& val){this->m_norms[5].set(i,j,val); }); 
      total_integrals += 6 ; 
    }
  }
  t_integral.stop();
  INFO( "Total integrals: " << total_integrals << " evaluate: " << t_eval << " integrate: " << t_integral );
  m_integralDispatch.flush();
}

double PolarisedAmplitude::prob( const AmpGen::Event& evt ) const {
  return m_weight * prob_unnormalised(evt) / m_norm;
}

void   PolarisedAmplitude::debug(const Event& evt ){
  for( auto& me : m_matrixElements )
  {
    INFO( me.decayTree->uniqueString() << " " << me.addressData 
        << " " << evt.getCache( me.addressData +0)
        << " " << evt.getCache( me.addressData +1 )
        << " " << evt.getCache( me.addressData +2)
        << " " << evt.getCache( me.addressData +3) );
  }
}

void PolarisedAmplitude::generateSourceCode( const std::string& fname, const double& normalisation, bool add_mt )
{
  INFO("Generating sourceCode -> " << fname );
  std::ofstream stream( fname );
  auto dim = m_eventType.dim(); 
  size_t size = dim.first * dim.second; 
  stream << "#include <complex>" << std::endl; 
  stream << "#include <vector>" << std::endl; 
  stream << "#include <array>" << std::endl; 
  Expression event = Parameter("x0",0,true,0); 
  std::vector<Array> perMatrixElement; 
  std::vector<Expression> expressions(size);
  for( auto& p : m_matrixElements ){
    p.pdf.prepare();
    p.pdf.to_stream( stream );
    p.pdf.compileWithParameters( stream );
    Array z( make_cse( Function( p.pdf.name() + "_wParams", {event} ) ), size );
    INFO( p.decayTree->uniqueString() << " coupling = " << p.coupling() );
    for( unsigned int j = 0 ; j < size; ++j ){
      expressions[j] = expressions[j] + p.coupling() * z[j];
    }
  }
  Tensor T_matrix( expressions, {dim.first,dim.second} );
  T_matrix.st();
  auto amplitude        = probExpression( T_matrix, { Constant(double(m_productionParameters[0])), 
                                                      Constant(double(m_productionParameters[1])) , 
                                                      Constant(double(m_productionParameters[2])) } );

  auto amplitude_extPol = probExpression( T_matrix, { Parameter("x2",0,true), 
                                                      Parameter("x3",0,true), 
                                                      Parameter("x4",0,true) } ); 
  stream << CompiledExpression< double, 
                                const double*, 
                                const int&>( amplitude / normalisation, "FCN",{},{}, m_mps ) << std::endl ;

  stream << CompiledExpression< double, 
                                const double*, 
                                const int&, 
                                const double&, 
                                const double&, 
                                const double&>( amplitude_extPol / normalisation, "FCN_extPol",{},{},m_mps ) << std::endl;
  stream.close();
}

Expression PolarisedAmplitude::probExpression( const Tensor& T_matrix, const std::vector<Expression>& p  ) const 
{
  Tensor T_conj = T_matrix.conjugate();
  Tensor::Index a,b,c; 
  Tensor TT = T_matrix(b,c) * T_conj(a,c);
  complex_t j(0,1);

  size_t it = T_matrix.dims()[0]; 
  Tensor rho( std::vector<size_t>({it,it}));
  if( it == 2 ){
    rho[{0,0}] = 1 + p[2]; 
    rho[{1,1}] = 1 - p[2];
    rho[{1,0}] = p[0] + j*p[1];
    rho[{0,1}] = p[0] - j*p[1];
  }
  else if ( it == 1 ) 
  {
    rho[{0,0}] = 1;
  }
  return Real( Expression( rho(a,b) * TT(b,a)  ));  
}

std::vector<FitFraction> PolarisedAmplitude::fitFractions( const LinearErrorPropagator& prop )
{
  return std::vector<FitFraction>();
}

void PolarisedAmplitude::transferParameters()
{ 
  m_probExpression.prepare(); 
  for( auto& me : m_matrixElements ) me.pdf.prepare();
}

real_t PolarisedAmplitude::getValNoCache( const AmpGen::Event& evt )  {
  transferParameters();
  AmpGen::Event copy(evt);
  copy.resizeCache( size() );
  for( auto& me :  m_matrixElements ){
    auto values = me(copy);
    copy.setCache( values , me.addressData ); 
  }
  return m_probExpression( copy.getCachePtr() );
}
