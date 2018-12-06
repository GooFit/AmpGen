#include "AmpGen/PolarisedSum.h"

#include <memory.h>
#include <ext/alloc_traits.h>
#include <chrono>
#include <complex>
#include <iomanip>
#include <map>
#include <memory>
#include <ostream>
#include <ratio>
#include <utility>

#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Array.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ThreadPool.h"

using namespace AmpGen;

PolarisedSum::PolarisedSum( const EventType& type, AmpGen::MinuitParameterSet& mps, const std::string& prefix ) : 
  m_mps(&mps),
  m_eventType(type)
{
  bool debug           = NamedParameter<bool>("PolarisedSum::Debug", false );
  bool autoCompile     = NamedParameter<bool>("PolarisedSum::AutoCompile",true);
  std::string objCache = NamedParameter<std::string>("PolarisedSum::ObjectCache",""); 

  AmplitudeRules proto( mps ); 

  auto proto_amplitudes = proto.getMatchingRules( type, prefix , true );
  auto production_polarisations = polarisations( type.mother() ); 
  std::vector<std::vector<int>> all_polarisation_states ;
  for( auto& pol : production_polarisations ) all_polarisation_states.push_back( {pol} );

  for( unsigned int i = 0 ; i < type.size(); ++i ){
    all_polarisation_states = 
      polarisationOuterProduct( all_polarisation_states, polarisations( type[i] ) );
  }
  auto set_polarisation_state = []( auto& matrix_element, auto& polState ){
    auto fs = matrix_element.first.getFinalStateParticles();
    matrix_element.first.setPolarisationState( polState[0] );
    for( unsigned int i = 0 ; i < fs.size(); ++i ){
      fs[i]->setPolarisationState( polState[i+1] );
    }
  };
  for( auto& m : proto_amplitudes ) INFO( m.first.uniqueString() ); 
  
  for( auto& matrix_element : proto_amplitudes ){
    Tensor thisExpression( std::vector<size_t>({all_polarisation_states.size()}) );
    int i = 0 ;
    DebugSymbols syms;  
    for( auto& polState : all_polarisation_states ){
      set_polarisation_state( matrix_element, polState );
      thisExpression[ i++ ] = make_cse( matrix_element.first.getExpression(&syms) ); 
    }
    CompiledExpression< std::vector<complex_t> , const real_t*, const real_t* > expression( 
        TensorExpression( thisExpression), 
        "p"+std::to_string( FNV1a_hash( matrix_element.first.uniqueString() ) ) , 
        type.getEventFormat(),debug ? syms : DebugSymbols() ,&mps ); 
    m_matrixElements.emplace_back( std::make_shared<Particle>(matrix_element.first), matrix_element.second, expression );
  } 
  for( auto& polState : all_polarisation_states){
    for( auto& matrix_element : proto_amplitudes ){
      set_polarisation_state( matrix_element, polState );
    }
  }
  m_polStates = all_polarisation_states; 
  if(autoCompile){
    ThreadPool tp(8);
    for( auto& thing : m_matrixElements ) 
      tp.enqueue( [&]{ CompilerWrapper().compile( thing.pdf, objCache ) ;}  );
  }
  if( mps.find("Px") == nullptr ){
    WARNING("Polarisation parameters not defined, defaulting to (0,0,0)");
  }
  m_pVector = { mps.addOrGet("Px",2,0,0) , mps.addOrGet("Py",2,0,0) , mps.addOrGet("Pz",2,0,0) };
  for( size_t i = 0 ; i < 6 ; ++i ) m_norms[i].resize( m_matrixElements.size(), m_matrixElements.size() );
}

std::vector<int> PolarisedSum::polarisations( const std::string& name ) const {
  auto props = *ParticlePropertiesList::get( name );
  std::vector<int> rt( props.twoSpin() + 1 );

  if( props.isFermion() ) return {1,-1};    
  if( props.twoSpin() == 0 ) return {0};
  if( props.twoSpin() == 2 ) 
    return (name == "gamma0") ? std::vector<int>({-1,1}) : std::vector<int>({-1,0,1}) ;
  else return {0};
}

std::vector<std::vector<int>> PolarisedSum::polarisationOuterProduct( const std::vector< std::vector< int > >& A, const std::vector< int >& B ) const {
  std::vector<std::vector<int>> rt; 
  for( auto& iA : A ){
    for( auto& iB : B ){
      rt.push_back( iA );
      rt.rbegin()->push_back( iB );
    }
  }
  return rt;
}

std::vector<TransitionMatrix<std::vector<complex_t>>> PolarisedSum::matrixElements() const
{
  return m_matrixElements;  
}

size_t count_zeros( std::array<Bilinears,6>& norms, const size_t& N )
{
  int nZeros = 0; 
  for( size_t i = 0 ; i < N; ++i ){
    for( size_t j = 0 ; j < N; ++j ){
      for( size_t k = 0 ; k < 6 ; ++k ) {
        double re = std::real( norms[k].get(i,j) );
        double im = std::imag( norms[k].get(i,j) );
        if( re < 1e-8 && im < 1e-8 ){ 
          nZeros++;
          norms[k].setZero(i,j);
        }
      }
    }
  }
  return nZeros; 
}

struct pclock {
  std::chrono::_V2::system_clock::time_point t_start; 
  std::chrono::_V2::system_clock::time_point t_end; 
  pclock() : t_start(   std::chrono::high_resolution_clock::now()) {}
  void stop(){ t_end = std::chrono::high_resolution_clock::now() ; }
  operator double() const { return std::chrono::duration<double, std::milli>( t_end - t_start ).count() ; } ; 
};

void   PolarisedSum::prepare()
{
  auto dim = m_eventType.dim();
  std::vector<size_t> changedPdfIndices; 
  pclock tEval; 
  for( size_t i = 0; i < m_matrixElements.size(); ++i ){
    pclock tMEval;
    auto& t = m_matrixElements[i];
    t.pdf.prepare();
    if( m_nCalls != 0 && !t.pdf.hasExternalsChanged() ) continue; 
    if( t.addressData == 999 ){
      DEBUG("Registering expression: " << t.decayTree->uniqueString() << " = " << t.addressData );
      t.addressData = m_events->registerExpression( t.pdf , dim.first * dim.second );
    }
    m_events->updateCache( t.pdf, t.addressData );
    tMEval.stop();
    t.pdf.resetExternals();
    changedPdfIndices.push_back(i);
    DEBUG("Updated cache for PDF: " << std::setw(55) << t.decayTree->uniqueString() << std::setw(3) << t.addressData << " " << "  t = " << tMEval << " ms" );
  }
  if( !m_probExpression.isLinked() ) build_probunnormalised();
  m_probExpression.prepare(); 
  m_weight = m_weightParam == nullptr ? 1 : m_weightParam->mean();
  tEval.stop(); 
  pclock tIntegral; 
  if( m_integrator.isReady() )
  {
    if( changedPdfIndices.size() != 0 ) calculateNorms(changedPdfIndices);
    std::array< complex_t, 6 > acc;
    for( size_t i = 0 ; i < m_matrixElements.size(); ++i )
    {
      for( size_t j = 0 ; j < m_matrixElements.size(); ++j )
      {
        complex_t c = m_matrixElements[i].coupling() * std::conj( m_matrixElements[j].coupling() );
        for( size_t k = 0 ; k < 6 ; ++k ) acc[k] += c * m_norms[k].get(i,j);
      }
    }
    double px = m_pVector[0];
    double py = m_pVector[1];
    double pz = m_pVector[2];
    complex_t z(px,-py);
    complex_t total =       (1+pz)*(acc[0]+acc[1]) +         (1-pz)*(acc[2]+acc[3]) 
                    + std::conj(z)*(acc[4]+acc[5]) + z * ( std::conj(acc[4]+acc[5]) );
    m_norm = std::real(total);
    //if( m_nCalls == 0 || changedPdfIndices.size() != 0 ) debug_norm();
  }
  tIntegral.stop();
  if( changedPdfIndices.size() != 0  ) 
    INFO("Time to evaluate = " << tEval << " ms; norm = " << tIntegral << " ms;  [pdfs = " << changedPdfIndices.size() << " zeros= " << count_zeros(m_norms, m_matrixElements.size() ) << "]" );
  m_nCalls++;
}

void PolarisedSum::debug_norm()
{
  double norm = 0;
  for( auto& evt : m_integrator.events() ){
    norm += evt.weight() * prob_unnormalised(evt) / evt.genPdf(); 
    if( prob_unnormalised(evt) > 1e6 ){
      WARNING("Event has very high weight: " << prob_unnormalised(evt) );
      evt.print();
    }
  }
  INFO("Average = " << m_norm << " " << norm /  m_integrator.sampleNorm()  );
  auto evt = m_integrator.events()[0];
  evt.print();
  INFO("Check one event: " << prob_unnormalised(evt) << " " << getValNoCache(evt) );
  int nNorms = 0 ; 
  int nZeros = 0; 
  for( size_t i = 0 ; i < m_matrixElements.size(); ++i ){
    for( size_t j = 0 ; j < m_matrixElements.size(); ++j ){
      std::string norm_string = "[";
      for( size_t k = 0 ; k < 6 ; ++k ) {
        double re = std::real( m_norms[k].get(i,j) );
        double im = std::imag( m_norms[k].get(i,j) );
        norm_string += "(" + std::to_string(re) +", " + std::to_string(im) + ") , ";
        if( re < 1e-8 && im < 1e-8 ) nZeros++;
      }
      norm_string += "]";
      nNorms += 6; 
      //INFO( m_matrixElements[i].decayTree->uniqueString() << " " << m_matrixElements[j].decayTree->uniqueString() << " " << norm_string );
    }
  }
  INFO( "nNorms = " << nNorms << " nZeros = " <<nZeros );
}

void   PolarisedSum::setEvents( AmpGen::EventList& events )
{ 
  reset();
  m_events = &events;
}

void   PolarisedSum::setMC( AmpGen::EventList& events )
{
  m_nCalls = 0;
  m_integrator = Integrator<18>(&events);
}

size_t PolarisedSum::size() const 
{ 
  auto dim = m_eventType.dim() ; 
  return dim.first * dim.second * m_matrixElements.size(); 
}

void   PolarisedSum::reset( const bool& flag ){ m_nCalls = 0 ; }

void   PolarisedSum::build_probunnormalised()
{
  Expression prob  = probExpression( transitionMatrix() , { Parameter("Px"), Parameter("Py"), Parameter("Pz") } );
  m_probExpression = CompiledExpression<real_t, const real_t*, const complex_t*>( 
      prob, "prob_unnormalised", std::map<std::string, size_t>(), {}, m_mps );
  CompilerWrapper().compile( m_probExpression );
} 
  
Tensor PolarisedSum::transitionMatrix()
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

complex_t PolarisedSum::TE( const Event& event, const size_t& x, const size_t& y )
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


double PolarisedSum::prob_unnormalised( const AmpGen::Event& evt ) const
{
  return m_probExpression( evt.getCachePtr(0) );
}

double PolarisedSum::norm() const 
{
  return m_norm;
}

void PolarisedSum::calculateNorms(const std::vector<size_t>& changedPdfIndices )
{
  size_t size_of = size() / m_matrixElements.size();
  std::vector<size_t> cacheIndex; 
  for( auto& i : changedPdfIndices ) m_integrator.prepareExpression( m_matrixElements[i].pdf , size_of );
  for( auto& m : m_matrixElements )  cacheIndex.push_back( m_integrator.events().getCacheIndex( m.pdf ) );
  
  for( auto& i : changedPdfIndices ){
  //for( auto i = 0 ; i < m_matrixElements.size(); ++i ){
    auto ai = cacheIndex[i]; 
    for(size_t j = 0; j < m_matrixElements.size(); ++j ){
      auto aj = cacheIndex[j];
      m_integrator.queueIntegral(ai+0, aj+0, i, j, &(m_norms[0]) );
      m_integrator.queueIntegral(ai+1, aj+1, i, j, &(m_norms[1]) );
      m_integrator.queueIntegral(ai+2, aj+2, i, j, &(m_norms[2]) );
      m_integrator.queueIntegral(ai+3, aj+3, i, j, &(m_norms[3]) ); 
      m_integrator.queueIntegral(ai+0, aj+2, i, j, &(m_norms[4]) , false );
      m_integrator.queueIntegral(ai+1, aj+3, i, j, &(m_norms[5]) , false );
      if( i != j )
      {
        m_integrator.queueIntegral(aj+0, ai+2, j, i, &(m_norms[4]) , false );
        m_integrator.queueIntegral(aj+1, ai+3, j, i, &(m_norms[5]) , false );
      }
    }
  }
  /*
  std::array< Bilinears,6 > copied_normalisations;
  for( int i = 0 ; i< 6 ; ++i ) 
    copied_normalisations[i].resize( m_matrixElements.size(), m_matrixElements.size() ); 
  
  for( auto i = 0 ; i < m_matrixElements.size(); ++i ){
    auto ai = cacheIndex[i]; 
    for(size_t j = 0; j < m_matrixElements.size(); ++j ){
      auto aj = cacheIndex[j];
      m_integrator.queueIntegral(ai+0, aj+0, i, j, &(copied_normalisations[0]) );
      m_integrator.queueIntegral(ai+1, aj+1, i, j, &(copied_normalisations[1]) );
      m_integrator.queueIntegral(ai+2, aj+2, i, j, &(copied_normalisations[2]) );
      m_integrator.queueIntegral(ai+3, aj+3, i, j, &(copied_normalisations[3]) ); 
      m_integrator.queueIntegral(ai+0, aj+2, i, j, &(copied_normalisations[4]) , false );
      m_integrator.queueIntegral(ai+1, aj+3, i, j, &(copied_normalisations[5]) , false );
    }
  }
  auto diff = [](const auto& l, const auto& r){
    return std::abs( l - r );
  };

  for( size_t i = 0 ; i < 6; ++i ){
   for( size_t j = 0 ; j < m_matrixElements.size(); ++j ){
     for( size_t k = 0 ; k < m_matrixElements.size(); ++k ){
       if( diff( m_norms[i](j,k) , copied_normalisations[i](j,k) ) > 1e-10 ) 
         INFO( "ERROR: normalisations( " << i  << ", " << j << ", " << k << " do not match !!"
             << m_norms[i](j,k) << " " << copied_normalisations[i](j,k)  <<
             diff( m_norms[i](j,k) , copied_normalisations[i](j,k))
             );
     }
   }
  }
  */
  m_integrator.flush();
  for( size_t i = 0 ; i < 6; ++i ) m_norms[i].resetCalculateFlags();
}

double PolarisedSum::prob( const AmpGen::Event& evt ) const {
  return m_weight * prob_unnormalised(evt) / m_norm;
}

void   PolarisedSum::debug(const Event& evt ){
  for( auto& me : m_matrixElements )
  {
    INFO( me.decayTree->uniqueString() << " " << me.addressData 
        << " " << evt.getCache(me.addressData+0)
        << " " << evt.getCache(me.addressData+1)
        << " " << evt.getCache(me.addressData+2)
        << " " << evt.getCache(me.addressData+3) );
  }
}

void PolarisedSum::generateSourceCode( const std::string& fname, const double& normalisation, bool add_mt )
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
  auto amplitude        = probExpression( T_matrix, { Constant(double(m_pVector[0])), 
                                                      Constant(double(m_pVector[1])) , 
                                                      Constant(double(m_pVector[2])) } );

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

Expression PolarisedSum::probExpression( const Tensor& T_matrix, const std::vector<Expression>& p  ) const 
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

std::vector<FitFraction> PolarisedSum::fitFractions( const LinearErrorPropagator& prop )
{
  return std::vector<FitFraction>();
}

void PolarisedSum::transferParameters()
{ 
  m_probExpression.prepare(); 
  for( auto& me : m_matrixElements ) me.pdf.prepare();
}

real_t PolarisedSum::getValNoCache( const AmpGen::Event& evt )  {
  transferParameters();
  AmpGen::Event copy(evt);
  copy.resizeCache( size() );
  for( auto& me :  m_matrixElements ){
    auto values = me(copy);
    copy.setCache( values , me.addressData ); 
  }
  return m_probExpression( copy.getCachePtr() );
}
