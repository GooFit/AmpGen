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
using namespace std::complex_literals; 

PolarisedSum::PolarisedSum( const EventType& type, 
                            MinuitParameterSet& mps, 
                            const std::string& prefix ) : 
  m_mps(&mps),
  m_eventType(type),
  m_prefix(prefix)
{
  m_debug              = NamedParameter<bool>(       "PolarisedSum::Debug"      ,false );
  std::string objCache = NamedParameter<std::string>("PolarisedSum::ObjectCache",""    );
  m_verbosity          = NamedParameter<bool>(       "PolarisedSum::Verbosity"  ,0     );
  m_rules              = AmplitudeRules(mps);
  auto protoAmps       = m_rules.getMatchingRules( type, prefix);
  auto prodPols        = polarisations( type.mother() ); 
  for(auto& pol : prodPols ) m_polStates.push_back({pol}); 
  for(size_t i = 0 ; i < type.size(); ++i ) m_polStates = polarisationOuterProduct(m_polStates, polarisations( type[i] ) );
  auto set_polarisation_state = []( auto& matrix_element, auto& polState ){
    auto fs = matrix_element.first.getFinalStateParticles();
    matrix_element.first.setPolarisationState(polState[0]);
    for(size_t i = 0 ; i < fs.size(); ++i ) fs[i]->setPolarisationState( polState[i+1] );
  };
  for( auto& m : protoAmps ) INFO( m.first.uniqueString() );  
  m_matrixElements.resize( protoAmps.size() );

  ThreadPool tp(8);
  for(unsigned i = 0; i < m_matrixElements.size(); ++i)
  {
    tp.enqueue( [i,this,&mps,&protoAmps, &set_polarisation_state]{
    Tensor thisExpression( Tensor::dim(this->m_polStates.size()) );
    DebugSymbols syms;  
    for( unsigned j=0;j<m_polStates.size();++j ){
      set_polarisation_state( protoAmps[i], m_polStates[j] );
      thisExpression[j] = make_cse( protoAmps[i].first.getExpression(&syms) ); 
    }
    m_matrixElements[i] = TransitionMatrix<std::vector<complex_t>>( 
        protoAmps[i].first,
        protoAmps[i].second, 
        CompiledExpression< std::vector<complex_t>,const real_t*, const real_t*>(
        TensorExpression(thisExpression), 
        protoAmps[i].first.decayDescriptor(),
        this->m_eventType.getEventFormat(), this->m_debug ? syms : DebugSymbols() ,&mps ) ); 
      CompilerWrapper().compile( m_matrixElements[i].amp );
    });
  } 
  
  auto   d = m_eventType.dim();
  auto   p = [&mps](const std::string& name){ return mps.addOrGet(name, Flag::Fix, 0, 0); };
  if( d.first == 1 )      m_pVector = {};
  else if( d.first == 2 ) m_pVector = {p("Px"), p("Py"), p("Pz")};
  else if( d.first == 3 ) m_pVector = {p("Px"), p("Py"), p("Pz"), p("Tyy"), p("Tzz"), p("Txy"), p("Txz"), p("Tyz")};
  size_t normSize = d.second * d.first * d.first;
  for(size_t i=0; i < normSize; ++i){
    m_norms.emplace_back( m_matrixElements.size(), m_matrixElements.size() );
  }  
}

std::vector<int> PolarisedSum::polarisations( const std::string& name ) const 
{
  auto props = *ParticlePropertiesList::get( name );
  if( props.twoSpin() == 0 ) return {0};
  if( props.twoSpin() == 1 ) return {1,-1};
  if( name == "gamma0" && props.twoSpin() == 2 ) return {1,-1};
  if( name != "gamma0" && props.twoSpin() == 2 ) return {1,0,-1};
  
  else { 
    WARNING("Particle with spin: " << props.twoSpin() << "/2" << " not implemented in initial/final state");
    return {0};
  }
}

std::vector<std::vector<int>> PolarisedSum::polarisationOuterProduct(const std::vector<std::vector<int>>& A, const std::vector<int>& B ) const 
{
  std::vector<std::vector<int>> rt; 
  for( auto& iA : A ){
    for( auto& iB : B ){
      rt.push_back(iA);
      rt.rbegin()->push_back(iB);
    }
  }
  return rt;
}

std::vector<complex_t> densityMatrix(const size_t& dim, const std::vector<MinuitProxy>& pv )
{
  if( dim == 1 ) return {1.};
  double px = pv[0];
  double py = pv[1];
  double pz = pv[2];
  if( dim == 2 ) return {1+pz    ,   px+1i*py,
                         px-1i*py,   1-pz     };
  if( dim == 3 ){
    double Tyy = pv[3];
    double Tzz = pv[4];
    double Txy = pv[5];
    double Txz = pv[6];
    double Tyz = pv[7];
    return {1 + 1.5*pz + sqrt(1.5)*Tzz                   , sqrt(0.375)*(px+1i*py) + sqrt(3.)*(Txz+1i*Tyz), -sqrt(1.5)*( Tzz + 2.*Tyy - 2.*1i*Txy), 
          sqrt(0.375)*(px-1i*py) + sqrt(3)*(Txz-1i*Tyz), 1 - sqrt(6.)*Tzz                              , sqrt(0.375)*(px+1i*py) - sqrt(3.)*(Txz+1i*Tyz) ,
         -sqrt(1.500)*( Tzz + 2.*Tyy + 2.*1i*Txy)      , sqrt(0.375)*(px-1i*py) - sqrt(3)*(Txz-1i*Tyz) , 1. - 1.5*pz + sqrt(1.5)*Tzz }; 
  }
  ERROR("Density matrices not implemented for state with size="<<dim);
  return {1.};
}

std::vector<Expression> convertProxies(const std::vector<MinuitProxy>& proxyVector, const std::function<Expression(const MinuitProxy&)>& transform)
{
  std::vector<Expression> rt;
  std::transform(proxyVector.begin(), proxyVector.end(), std::back_inserter(rt), transform );
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
  std::vector<bool> hasChanged( m_matrixElements.size(), false); 
  size_t nChanges = 0; 
  ProfileClock tEval; 
  size_t size_of = size() / m_matrixElements.size();
  if( m_events != nullptr ) m_events->reserveCache( size() );
  if( m_integrator.isReady() ) m_integrator.events().reserveCache( size() );
  for( size_t i = 0; i < m_matrixElements.size(); ++i ){
    ProfileClock tMEval;
    auto& t = m_matrixElements[i];
    if( m_nCalls != 0 && !t.amp.hasExternalsChanged() ) continue; 
    if( t.addressData == 999 ) t.addressData = m_events->registerExpression(t.amp, dim.first * dim.second );
    m_events->updateCache(t.amp, t.addressData);
    m_integrator.prepareExpression(t.amp, size_of);
    tMEval.stop();
    t.amp.resetExternals();
    hasChanged[i] = true; 
    nChanges++;
    if( m_nCalls == 0 && m_integrator.isReady() ) m_integIndex.push_back( m_integrator.events().getCacheIndex( t.amp ) );
  }
  if( !m_probExpression.isLinked() ) build_probunnormalised();
  m_weight = m_weightParam == nullptr ? 1 : m_weightParam->mean();
  tEval.stop(); 
  ProfileClock tIntegral;  
  m_rho = densityMatrix(m_eventType.dim().first, m_pVector);
  if( m_integrator.isReady() )
  {
    if(nChanges != 0) calculateNorms(hasChanged);
    complex_t z = 0;
    for(size_t i = 0; i < m_matrixElements.size(); ++i){
      for(size_t j = i; j < m_matrixElements.size(); ++j){
         z += ((i==j) ? 1. : 2. ) * m_matrixElements[i].coupling()*std::conj(m_matrixElements[j].coupling())*norm(i,j);
      }
    }
    m_norm = std::real(z); 
    if(m_nCalls % 10000 == 0 && m_prefix == "") debug_norm();
  }
  tIntegral.stop();
  if(m_verbosity && nChanges != 0)
    INFO("Time to evaluate = " << tEval << " ms; " << "norm = "  << tIntegral << " ms; " << "pdfs = "  << nChanges);
  m_nCalls++;
}

void PolarisedSum::debug_norm()
{
  double norm_slow = 0;
  for( auto& evt : m_integrator.events() ) 
    norm_slow += evt.weight() * prob_unnormalised(evt) / evt.genPdf();
  auto evt = m_integrator.events()[0];
  INFO("Event[0]: " << prob_unnormalised(evt) << " " << getValNoCache(evt) );
  INFO("Norm    : " << std::setprecision(10) 
      << "bilinears=" << m_norm 
      << "; exact="     << norm_slow /  m_integrator.sampleNorm()
      << "; d = " << m_norm - norm_slow / m_integrator.sampleNorm()  
      << "; sample=" << m_integrator.sampleNorm() );
}

void   PolarisedSum::setEvents( EventList& events )
{ 
  reset();
  m_events = &events;
}

void   PolarisedSum::setMC( EventList& events )
{
  m_nCalls = 0;
  m_integrator = integrator(&events);
}

size_t PolarisedSum::size() const 
{ 
  auto dim = m_eventType.dim() ; 
  return dim.first * dim.second * m_matrixElements.size(); 
}

void   PolarisedSum::reset( const bool& flag ){ m_nCalls = 0 ; }

void   PolarisedSum::build_probunnormalised()
{
  DebugSymbols db; 
  auto prob = probExpression(transitionMatrix(), convertProxies(m_pVector,[](auto& p){ return Parameter(p->name());} ), m_debug ? &db : nullptr);
  m_probExpression = CompiledExpression<real_t, const real_t*, const complex_t*>(prob, "prob_unnormalised", std::map<std::string, size_t>(), db, m_mps);
  CompilerWrapper().compile(m_probExpression);
  m_probExpression.prepare();
} 
  
Tensor PolarisedSum::transitionMatrix()
{ 
  auto dim = m_eventType.dim();
  auto size = dim.first * dim.second ; 
  std::vector<Expression> expressions(size, 0);
  for( auto& me : m_matrixElements ){
    auto coupling   = me.coupling.to_expression() ;
    auto cacheIndex = m_events->getCacheIndex(me.amp); 
    for( size_t i = 0 ; i < size ; ++i ){
      expressions[i] = expressions[i] + coupling * Parameter( "x1["+std::to_string(cacheIndex+i)+"]",0,true); 
    }
  }
  Tensor T_matrix(expressions, {dim.first, dim.second});
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

complex_t PolarisedSum::norm(const size_t& i, const size_t& j, AmpGen::PolarisedSum::integrator* integ)
{
  auto   ai = m_integIndex[i];
  auto   aj = m_integIndex[j];
  complex_t total = 0;
  auto   dim = m_eventType.dim();
  auto    s1 = dim.first;
  auto    s2 = dim.second;
  for(size_t x = 0 ; x < m_norms.size(); ++x){
    auto f         = x % s2;
    auto psiIndex  = (x-f) / s2;
    auto m2        = psiIndex % s1;
    auto m1        = (psiIndex-m2)/s1;
    total          += m_rho[psiIndex] * m_norms[x].get(i, j, integ, ai+m1*s2+f, aj+m2*s2+f);
  }
  return total; 
}

void PolarisedSum::calculateNorms(const std::vector<bool>& hasChanged)
{
  for( size_t i = 0 ; i < m_matrixElements.size(); ++i ){
    for( size_t j = i; j < m_matrixElements.size(); ++j ){
      if( hasChanged[i] || hasChanged[j] ) norm(i, j, &m_integrator);
    }
  }
  m_integrator.flush();
}

double PolarisedSum::prob(const Event& evt) const 
{
  return m_weight * prob_unnormalised(evt) / m_norm;
}

void PolarisedSum::debug(const Event& evt)
{
  auto dim    = m_eventType.dim(); 
  size_t tsize = dim.first * dim.second;   
  for(auto& me : m_matrixElements)
  {
    std::vector<complex_t> this_cache(0,tsize);
    for(size_t i = 0 ; i < tsize; ++i ) this_cache.emplace_back( evt.getCache(me.addressData+i) );
    INFO( me.decayDescriptor() << " " << vectorToString( this_cache, " ") );
  }
  INFO("P(x) = " << getValNoCache(evt) );
  INFO("Prod = [" << vectorToString(m_pVector , ", ") <<"]");
  if( m_debug )
  {
    transferParameters();
    Event copy(evt);
    copy.resizeCache( size() );
    for(auto& me : m_matrixElements){
      auto values = me(copy);
      copy.setCache( values , me.addressData );
    }
    m_probExpression.debug( copy.getCachePtr() );
  }
}

void PolarisedSum::generateSourceCode(const std::string& fname, const double& normalisation, bool add_mt)
{
  INFO("Generating sourceCode -> " << fname );
  std::ofstream stream( fname );
  auto dim = m_eventType.dim(); 
  size_t size = dim.first * dim.second; 
  CompilerWrapper().preamble( stream );
  Expression event = Parameter("x0",0,true); 
  std::vector<Expression> expressions(size);
  for( auto& p : m_matrixElements ){
    p.amp.prepare();
    p.amp.to_stream( stream );
    p.amp.compileWithParameters( stream );
    Array z( make_cse( Function( programatic_name( p.amp.name()) + "_wParams", {event} ) ), size );
    INFO( p.decayDescriptor() << " coupling = " << p.coupling() );
    for( unsigned int j = 0 ; j < size; ++j ){
      expressions[j] = expressions[j] + p.coupling() * z[j];
    }
  }
  Tensor T_matrix( expressions, {dim.first,dim.second} );
  T_matrix.st();
  auto amp        = probExpression(T_matrix, convertProxies(m_pVector, [](auto& proxy) -> Expression{ return double(proxy);} )); 
  auto amp_extPol = probExpression(T_matrix, {Parameter("x2",0,true), Parameter("x3",0,true), Parameter("x4",0,true)}); 
  stream << CompiledExpression<double, 
                               const double*, 
                               const int&>( amp / normalisation, "FCN",{},{}, m_mps ) << std::endl ;

  stream << CompiledExpression<double, 
                               const double*, 
                               const int&, 
                               const double&, 
                               const double&, 
                               const double&>( amp_extPol / normalisation, "FCN_extPol",{},{},m_mps ) << std::endl;
  stream.close();
}

Expression PolarisedSum::probExpression(const Tensor& T_matrix, const std::vector<Expression>& p, DebugSymbols* db) const 
{
  Tensor T_conj = T_matrix.conjugate();
  Tensor::Index a,b,c; 
  Tensor TT = T_matrix(a,b) * T_conj(c,b);
  size_t it = T_matrix.dims()[0]; 
  Tensor rho = Identity(it);
  if(it == 2) rho = rho + Sigma[0] * p[0] + Sigma[1] * p[1] + Sigma[2]*p[2];  
  if(it == 3)
  { 
    auto px  = p[0];
    auto py  = p[1];
    auto pz  = p[2];
    auto Tyy = p[3];
    auto Tzz = p[4];
    auto Txy = p[5];
    auto Txz = p[6];
    auto Tyz = p[7];
    rho(0,0) = 1 + 1.5 * pz + sqrt(1.5)*Tzz;
    rho(1,0) = sqrt(0.375)*(px+1i*py) + sqrt(3.)*(Txz+1i*Tyz);
    rho(2,0) = -sqrt(1.5)*( Tzz + 2.*Tyy - 2.*1i*Txy); 
    rho(0,1) = sqrt(0.375)*(px-1i*py) + sqrt(3.)*(Txz-1i*Tyz);
    rho(1,1) = 1 - sqrt(6.)*Tzz;
    rho(2,1) = sqrt(0.375)*(px+1i*py) - sqrt(3.)*(Txz+1i*Tyz);
    rho(0,2) = -sqrt(1.5)*( Tzz + 2.*Tyy + 2.*1i*Txy); 
    rho(1,2) = sqrt(0.375)*(px-1i*py) - sqrt(3)*(Txz-1i*Tyz);
    rho(2,2) = 1. - 1.5*pz + sqrt(1.5)*Tzz;
  }
  ADD_DEBUG_TENSOR(T_matrix, db);
  ADD_DEBUG_TENSOR(rho, db);
  ADD_DEBUG_TENSOR(TT , db);
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
    me.amp.prepare();
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

void PolarisedSum::setWeight( MinuitParameter* param ){ m_weightParam = param ; } 
double PolarisedSum::getWeight() const { return m_weightParam == nullptr ? 1.0 : m_weightParam->mean() ; }
