#include "AmpGen/PolarisedSum.h"

#include <memory.h>
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
#include "AmpGen/ProfileClock.h"
#include "AmpGen/DiracMatrices.h"
#include "AmpGen/Simplify.h"

using namespace AmpGen;

PolarisedSum::PolarisedSum( const EventType& type, 
                            MinuitParameterSet& mps, 
                            const std::string& prefix ) : 
  m_mps(&mps),
  m_eventType(type),
  m_prefix(prefix)
{
  bool debug           = NamedParameter<bool>(       "PolarisedSum::Debug"      ,false );
  bool autoCompile     = NamedParameter<bool>(       "PolarisedSum::AutoCompile",true  );
  std::string objCache = NamedParameter<std::string>("PolarisedSum::ObjectCache",""    );
  m_verbosity          = NamedParameter<bool>(       "PolarisedSum::Verbosity"  ,0     );
  m_rules = AmplitudeRules(mps);

  auto proto_amplitudes = m_rules.getMatchingRules( type, prefix);
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
        matrix_element.first.decayDescriptor(),
        type.getEventFormat(),debug ? syms : DebugSymbols() ,&mps ); 
    m_matrixElements.emplace_back(matrix_element.first, matrix_element.second, expression );
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
  if( mps.find("Px") == nullptr ) WARNING("Polarisation parameters not defined, defaulting to (0,0,0)");
  m_pVector = {mps.addOrGet("Px",2,0,0), mps.addOrGet("Py",2,0,0), mps.addOrGet("Pz",2,0,0)};
  
  auto   d = m_eventType.dim();
  size_t normSize = d.second * d.first * ( d.first + 1)/2;
  for(size_t i=0; i < normSize; ++i) m_norms.emplace_back( m_matrixElements.size(), m_matrixElements.size() );
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

std::vector<std::vector<int>> PolarisedSum::polarisationOuterProduct( const std::vector<std::vector<int>>& A, const std::vector<int>& B ) const {
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

void   PolarisedSum::prepare()
{
  DEBUG( "Preparing: " << m_prefix << " " << m_events << " ready = " << m_integrator.isReady() );
  transferParameters();
  auto dim = m_eventType.dim();
  std::vector<size_t> changedPdfIndices; 
  ProfileClock tEval; 
  size_t size_of = size() / m_matrixElements.size();
  for( size_t i = 0; i < m_matrixElements.size(); ++i ){
    ProfileClock tMEval;
    auto& t = m_matrixElements[i];
    if( m_nCalls != 0 && !t.pdf.hasExternalsChanged() ) continue; 
    if( t.addressData == 999 ){
    //  INFO("Registering expression: " << t.decayDescriptor() << " = " << t.addressData );
      t.addressData = m_events->registerExpression( t.pdf , dim.first * dim.second );
    }
    m_events->updateCache(t.pdf, t.addressData);
    m_integrator.prepareExpression(t.pdf, size_of);
    tMEval.stop();
    t.pdf.resetExternals();
    changedPdfIndices.push_back(i);
   // INFO("Updated cache for PDF: " << std::setw(55) << t.decayDescriptor() << std::setw(3) << t.addressData << " " << "  t = " << tMEval << " ms" );
  }
  if( !m_probExpression.isLinked() ) build_probunnormalised();
  m_weight = m_weightParam == nullptr ? 1 : m_weightParam->mean();
  tEval.stop(); 
  ProfileClock tIntegral; 
  if( m_integrator.isReady() )
  {
    DEBUG("size(transitionMatrix) = " << m_norms.size() << "; size(matrixElements) = " << m_matrixElements.size() );
    if( changedPdfIndices.size() != 0 ) calculateNorms(changedPdfIndices);
    std::vector<complex_t> acc(m_norms.size());
    for(size_t i = 0; i < m_matrixElements.size(); ++i)
    {
      for(size_t j = 0; j < m_matrixElements.size(); ++j)
      {
        complex_t c = m_matrixElements[i].coupling() * std::conj( m_matrixElements[j].coupling() );
        for(size_t k = 0; k < m_norms.size(); ++k) acc[k] += c * m_norms[k].get(i,j);
      }
    }
    m_norm = std::real(acc[0]+acc[1]+acc[2]+acc[3]) 
           + std::real(acc[0]+acc[1]-acc[2]-acc[3]) * m_pVector[2]
       + 2 * std::real((acc[4]+acc[5])*complex_t(m_pVector[0],m_pVector[1]));
//    if(m_nCalls == 0 && m_prefix == "") debug_norm();
  }
  tIntegral.stop();
  if(m_verbosity && changedPdfIndices.size() != 0)
    INFO("Time to evaluate = " << tEval                  << " ms; "
      << "norm = "  << tIntegral << " ms; "
      << "pdfs = "  << changedPdfIndices.size() );
  m_nCalls++;
}

void PolarisedSum::debug_norm()
{
  double norm_slow = 0;
  for( auto& evt : m_integrator.events() ) norm_slow += evt.weight() * prob_unnormalised(evt) / evt.genPdf();

  INFO("NORM = " << std::setprecision(10) << " " << m_norm << " " << norm_slow /  m_integrator.sampleNorm() << " sample norm = " << m_integrator.sampleNorm() );
  auto evt = m_integrator.events()[0];
  INFO("[" << m_prefix << "]  Check one event: " << std::setprecision(10) << " " << prob_unnormalised(evt) << " " << getValNoCache(evt) );
  int nNorms       = 0;
  int nZeros       = 0; 
  int nConjugate   = 0;
  int nNegatives   = 0; 
  int nReplicas    = 0; 
  int negConjugate = 0;
  auto s = m_matrixElements.size(); 
  std::vector<bool> flagged( 6 * s*s, false ); 
  auto identify = [&flagged,&s]( const size_t i2, const size_t j2, const size_t& k2, int& counter ){
    flagged[k2 + s * j2 + s*s*i2] = true;
    counter++;
  }; 
  complex_t t(0,0);
  for( size_t i = 0 ; i < s; ++i ){
    for( size_t j = 0 ; j < s; ++j ){
      complex_t c = m_matrixElements[i].coupling() * std::conj( m_matrixElements[j].coupling() );
      t += c * norm(i,j);
      for( size_t k = 0 ; k < 6 ; ++k ) {
        if( flagged[ k + s*j + s*s*i] ) continue; 
        double re = std::real( m_norms[k].get(i,j) );
        double im = std::imag( m_norms[k].get(i,j) );
        if( re < 1e-8 && im < 1e-8 ) identify(i,j,k,nZeros);
        for( size_t i2 = i; i2 < m_matrixElements.size(); ++i2 ){
          for( size_t j2 = j; j2 < m_matrixElements.size(); ++j2 ){
            for( size_t k2 = 0 ; k2 < 6 ; ++k2 ) {
              if( i2 == i && j2 == j && k2 == k ) continue;
              if( flagged[k2 + s*j2 + s*s*i2] ) continue; 
              double re2 = std::real( m_norms[k2].get(i2,j2) );
              double im2 = std::imag( m_norms[k2].get(i2,j2) );
              if( std::abs(re2 - re) < 1e-6 && std::abs(im2 - im) < 1e-6 ) identify(i2,j2,k2, nReplicas );
              if( std::abs(re2 - re) < 1e-6 && std::abs(im2 + im) < 1e-6 ) identify(i2,j2,k2, nConjugate );
              if( std::abs(re2 + re) < 1e-6 && std::abs(im2 + im) < 1e-6 ) identify(i2,j2,k2, nNegatives );
              if( std::abs(re2 + re) < 1e-6 && std::abs(im2 - im) < 1e-6 ) identify(i2,j2,k2, negConjugate );
            }  
          }
        }
      }
      nNorms += 6; 
    }
  }
  INFO( "nNorms = " << nNorms << " nZeros = " <<nZeros << " copies = " << nReplicas << " conj = " << nConjugate << " nNegatives = " << nNegatives << " negConj = " << negConjugate );
  INFO( std::real(t));
}

void   PolarisedSum::setEvents( EventList& events )
{ 
  reset();
  m_events = &events;
}

void   PolarisedSum::setMC( EventList& events )
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
  m_probExpression.prepare();
} 
  
Tensor PolarisedSum::transitionMatrix()
{ 
  auto dim = m_eventType.dim();
  auto size = dim.first * dim.second ; 
  std::vector<Expression> expressions( size, 0);
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

double PolarisedSum::prob_unnormalised( const Event& evt ) const
{
  return m_probExpression( evt.getCachePtr(0) );
}

double PolarisedSum::norm() const 
{
  return m_norm;
}

complex_t PolarisedSum::norm(const size_t& i, const size_t& j) const 
{
  double px = m_pVector[0];
  double py = m_pVector[1];
  double pz = m_pVector[2];
  return (1+pz)*(m_norms[0].get(i,j)+m_norms[1].get(i,j)) 
          + (1-pz)*(m_norms[2].get(i,j)+m_norms[3].get(i,j))
          + complex_t(px, py) * ( m_norms[4].get(i,j) + m_norms[5].get(i,j) ) 
          + complex_t(px,-py) * std::conj( m_norms[4].get(j,i) + m_norms[5].get(j,i) );
}

void PolarisedSum::calculateNorms(const std::vector<size_t>& changedPdfIndices)
{
  std::vector<size_t> cacheIndex; 
  for(auto& m : m_matrixElements ) cacheIndex.push_back(m_integrator.events().getCacheIndex(m.pdf));
  for(auto& i : changedPdfIndices){
    auto ai = cacheIndex[i]; 
    for(size_t j = 0; j < m_matrixElements.size(); ++j){
      auto aj = cacheIndex[j];
      m_integrator.queueIntegral(ai+0, aj+0, i, j, &(m_norms[0]));
      m_integrator.queueIntegral(ai+1, aj+1, i, j, &(m_norms[1]));
      m_integrator.queueIntegral(ai+2, aj+2, i, j, &(m_norms[2]));
      m_integrator.queueIntegral(ai+3, aj+3, i, j, &(m_norms[3])); 
      m_integrator.queueIntegral(ai+0, aj+2, i, j, &(m_norms[4]), false);
      m_integrator.queueIntegral(ai+1, aj+3, i, j, &(m_norms[5]), false);
      if( i != j )
      {
        m_integrator.queueIntegral(aj+0, ai+2, j, i, &(m_norms[4]), false);
        m_integrator.queueIntegral(aj+1, ai+3, j, i, &(m_norms[5]), false);
      }
    }
  }
  m_integrator.flush();
  for(auto& norm : m_norms) norm.resetCalculateFlags();
}

double PolarisedSum::prob(const Event& evt) const 
{
  return m_weight * prob_unnormalised(evt) / m_norm;
}

void PolarisedSum::debug(const Event& evt)
{
  auto dim = m_eventType.dim(); 
  size_t size = dim.first * dim.second; 
  
  for(auto& me : m_matrixElements)
  {
    std::vector<complex_t> this_cache(0,size);
    for( int i = 0 ; i < size; ++i ) this_cache.emplace_back( evt.getCache(me.addressData+i) );
    INFO( me.decayDescriptor() << " " << vectorToString( this_cache, " ") );
  }
}

void PolarisedSum::generateSourceCode(const std::string& fname, const double& normalisation, bool add_mt)
{
  INFO("Generating sourceCode -> " << fname );
  std::ofstream stream( fname );
  auto dim = m_eventType.dim(); 
  size_t size = dim.first * dim.second; 
  CompilerWrapper().preamble( stream );
  Expression event = Parameter("x0",0,true,0); 
  //std::vector<Array> perMatrixElement; 
  std::vector<Expression> expressions(size);
  for( auto& p : m_matrixElements ){
    p.pdf.prepare();
    p.pdf.to_stream( stream );
    p.pdf.compileWithParameters( stream );
    Array z( make_cse( Function( programatic_name( p.pdf.name()) + "_wParams", {event} ) ), size );

    INFO( p.decayDescriptor() << " coupling = " << p.coupling() );
    for( unsigned int j = 0 ; j < size; ++j ){
      expressions[j] = expressions[j] + p.coupling() * z[j];
    }
  }
  Tensor T_matrix( expressions, {dim.first,dim.second} );
  T_matrix.st();
  auto amplitude        = probExpression(T_matrix, {Constant(double(m_pVector[0])), 
                                                    Constant(double(m_pVector[1])), 
                                                    Constant(double(m_pVector[2]))});

  auto amplitude_extPol = probExpression(T_matrix, {Parameter("x2",0,true), 
                                                    Parameter("x3",0,true), 
                                                    Parameter("x4",0,true)}); 
  stream << CompiledExpression<double, 
                               const double*, 
                               const int&>( amplitude / normalisation, "FCN",{},{}, m_mps ) << std::endl ;

  stream << CompiledExpression<double, 
                               const double*, 
                               const int&, 
                               const double&, 
                               const double&, 
                               const double&>( amplitude_extPol / normalisation, "FCN_extPol",{},{},m_mps ) << std::endl;
  stream.close();
}

Expression PolarisedSum::probExpression(const Tensor& T_matrix, const std::vector<Expression>& p) const 
{
  Tensor T_conj = T_matrix.conjugate();
  Tensor::Index a,b,c; 
  Tensor TT = T_matrix(a,b) * T_conj(c,b);
  size_t it = T_matrix.dims()[0]; 
  Tensor rho = Identity(it);
  if(it == 2) rho = rho + Sigma[0] * p[0] + Sigma[1] * p[1] + Sigma[2]*p[2];
  Expression rt = rho(a,b) * TT(b,a);
  return Real(rt);  
}

std::vector<FitFraction> PolarisedSum::fitFractions(const LinearErrorPropagator& prop)
{
  bool recomputeIntegrals    = NamedParameter<bool>("PolarisedSum::RecomputeIntegrals", false );
  std::vector<FitFraction> outputFractions;
 
  for(auto& rule : m_rules.rules()) 
  {
    FitFractionCalculator<PolarisedSum> pCalc(this, findIndices(m_matrixElements, rule.first), recomputeIntegrals);
    INFO("Denom = [" << vectorToString( findIndices(m_matrixElements, rule.first ) , ", " ) << "]");
    for(auto& process : rule.second) 
    {
      if(process.head() == m_eventType.mother() && process.prefix() != m_prefix) continue;
      auto numeratorIndices   = processIndex(m_matrixElements, process.name());
      if(numeratorIndices.size() == 0 || numeratorIndices == pCalc.normSet ) continue; 
      INFO("Adding calculation: " << process.name() << " [" << vectorToString(numeratorIndices, ", ") << "]");    
      pCalc.emplace_back(process.name(), numeratorIndices);
    }
    if( pCalc.calculators.size() == 0 ) continue;  
    auto fractions = pCalc(rule.first, prop);
    for( auto& f : fractions ) outputFractions.emplace_back(f);
  }
  auto head_rules = m_rules.rulesForDecay(m_eventType.mother(), m_prefix); 
  
  FitFractionCalculator<PolarisedSum> iCalc(this, findIndices(m_matrixElements, m_eventType.mother()), recomputeIntegrals);
  for(size_t i = 0 ; i < head_rules.size(); ++i)
  {
    auto process_i = head_rules[i];
    auto num_i   = processIndex(m_matrixElements, process_i.name());
    if( num_i.size() == 0 || num_i == iCalc.normSet ) continue; 
    for( size_t j = i+1 ; j < head_rules.size(); ++j ){
      auto process_j = head_rules[j];
      auto num_j   = processIndex(m_matrixElements, process_j.name());
      if( num_j.size() == 0 || num_j == iCalc.normSet ) continue; 
      iCalc.emplace_back(process_i.name() + " " + process_j.name() , num_i, num_j);
    }
  }
  auto ifractions = iCalc(m_eventType.mother(), prop);
  for(auto& p : outputFractions) INFO(p);
  
  INFO("INTERFERENCE FRACTIONS");
  for( auto& f : ifractions ) INFO( FitFraction(f) );
  return outputFractions;
}

void PolarisedSum::transferParameters()
{ 
  if( m_probExpression.isLinked() ) m_probExpression.prepare();
  for(auto& me : m_matrixElements){
    me.coefficient = me.coupling();
    me.pdf.prepare();
  }
  for(auto& p  : m_pVector ) p.update();
}

real_t PolarisedSum::getValNoCache( const Event& evt )
{
  transferParameters();
  Event copy(evt);
  copy.resizeCache( size() );
  for(auto& me : m_matrixElements){
    auto values = me(copy);
    copy.setCache( values , me.addressData ); 
  }
  return m_probExpression( copy.getCachePtr() );
}
