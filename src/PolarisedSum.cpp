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
#include <thread>

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
#include "AmpGen/enum.h"
#include "AmpGen/simd/utils.h"

using namespace AmpGen;
using namespace std::complex_literals; 

#if DEBUGLEVEL == 1  
  ENABLE_DEBUG( PolarisedSum )
#endif

namespace AmpGen { make_enum(spaceType, spin, flavour) }

std::vector<Expression> convertProxies(const std::vector<MinuitProxy>& proxyVector, const std::function<Expression(const MinuitProxy&)>& transform)
{
  std::vector<Expression> rt;
  std::transform(proxyVector.begin(), proxyVector.end(), std::back_inserter(rt), transform );
  return rt; 
}

PolarisedSum::PolarisedSum(const EventType& type, 
                           MinuitParameterSet& mps, 
                           const std::vector<MinuitProxy>& pVector) 
  : m_mps       (&mps)
  , m_pVector   (pVector)
  , m_verbosity (NamedParameter<bool>("PolarisedSum::Verbosity", 0     ))
  , m_debug     (NamedParameter<bool>("PolarisedSum::Debug"    , false ))
  , m_eventType (type)
  , m_rules     (mps)
  , m_dim       (m_eventType.dim())
{
  std::string objCache = NamedParameter<std::string>("PolarisedSum::ObjectCache", ""    );
  spaceType stype      = NamedParameter<spaceType>(  "PolarisedSum::SpaceType"  , spaceType::spin);
  {
  ThreadPool tp(std::thread::hardware_concurrency() );
  if( stype == spaceType::spin )
  {
    auto prodPols        = polarisations(m_eventType.mother()); 
    std::vector<std::vector<int>> polStates; 
    for(const auto& pol : prodPols ) polStates.push_back({pol}); 
    for(unsigned i = 0 ; i != type.size(); ++i ) polStates = indexProduct(polStates, polarisations( type[i] ) ); 
    auto protoAmps       = m_rules.getMatchingRules(m_eventType);
    for(const auto& m : protoAmps ) INFO( m.first.uniqueString() ); 
    m_matrixElements.resize( protoAmps.size() );
    for(unsigned i = 0; i < m_matrixElements.size(); ++i)
    {
      auto [lp, lc] = protoAmps[i];
      auto & p = lp;
      auto & c = lc;
      PolarisedSum* ptr = this; 
      tp.enqueue( [i, p=lp, c=lc, polStates, &mps, ptr] () mutable {
        Tensor thisExpression( Tensor::dim(polStates.size()) );
        DebugSymbols syms;      
        for(unsigned j = 0; j != polStates.size(); ++j) 
          thisExpression[j] = make_cse( p.getExpression(j == 0 ? &syms: nullptr, polStates[j] ) );         
          //thisExpression[j] = make_cse( p.getExpression(&syms, polStates[j] ) );         
        ptr->m_matrixElements[i] = TransitionMatrix<void>( 
            p, c,
            CompiledExpression<void(complex_v*, const size_t&, const real_t*, const float_v*)>(
            TensorExpression(thisExpression), p.decayDescriptor(), &mps,
            ptr->m_eventType.getEventFormat(), ptr->m_debug ? syms : DebugSymbols() ) );
        
        CompilerWrapper().compile( ptr->m_matrixElements[i] );
        ptr->m_matrixElements[i].size = thisExpression.size();
      });
    }
  }
  if ( stype == spaceType::flavour )
  {
    m_dim = {2,1};
    auto r1 = m_rules.getMatchingRules(m_eventType, m_prefix);
    auto r2 = m_rules.getMatchingRules(m_eventType.conj(true), m_prefix);  
    m_matrixElements.resize( r1.size() + r2.size() );
    for(unsigned i = 0 ; i != m_matrixElements.size(); ++i)
    {
      tp.enqueue( [i, this, &r1, &r2] () mutable {
        Tensor thisExpression( Tensor::dim(2) );
        DebugSymbols syms;  
        auto& [p,coupling] = i < r1.size() ? r1[i] : r2[i-r1.size()];
        thisExpression[0] = i < r1.size() ? make_cse( p.getExpression(&syms) ) : 0;
        thisExpression[1] = i < r1.size() ? 0 : make_cse( p.getExpression(&syms) ); 
        this->m_matrixElements[i] = TransitionMatrix<void>( 
            p, coupling, 
            CompiledExpression<void(complex_v*, const size_t&, const real_t*, const float_v*)>(
            TensorExpression(thisExpression), p.decayDescriptor(), this->m_mps,
            this->m_eventType.getEventFormat(), this->m_debug ? syms : DebugSymbols() ) );
          CompilerWrapper().compile( m_matrixElements[i] );
        });
    }
  }
  }
  if( m_pVector.size() == 0 )
  {
    auto p = [this](const std::string& name){ return this->m_mps->addOrGet(name, Flag::Fix, 0, 0); };
    if( m_dim.first == 1 )      m_pVector = {};
    else if( m_dim.first == 2 ) m_pVector = {p("Px"), p("Py"), p("Pz")};
    else if( m_dim.first == 3 ) m_pVector = {p("Px"), p("Py"), p("Pz"), p("Tyy"), p("Tzz"), p("Txy"), p("Txz"), p("Tyz")};
  }
  for(size_t i=0; i < m_dim.second * m_dim.first * m_dim.first; ++i) m_norms.emplace_back( m_matrixElements.size(), m_matrixElements.size() ); 
  
  DebugSymbols db; 
  auto prob = probExpression(transitionMatrix(), convertProxies(m_pVector,[](auto& p){ return Parameter(p->name());} ), m_debug ? &db : nullptr);
  m_probExpression = make_expression<float_v, real_t, complex_v>( prob, "prob_unnormalised", m_mps, this->m_debug ? db : DebugSymbols() );
}

std::vector<int> PolarisedSum::polarisations( const std::string& name ) const 
{
  auto props = *ParticlePropertiesList::get( name );
  if( props.twoSpin() == 0 ) return {0};         // scalar 
  if( props.isPhoton() )     return {1,-1};      // photon 
  if( props.twoSpin() == 1 ) return {1,-1};      // fermion
  if( props.twoSpin() == 4 ) return {-2,1,0,1,2};// tensor
  if( props.twoSpin() == 2 ) return {1,0,-1};    // vector
  else { 
    WARNING("Particle with spin: " << props.twoSpin() << "/2" << " not implemented in initial/final state");
    return {0};
  }
}

std::vector<std::vector<int>> PolarisedSum::indexProduct(const std::vector<std::vector<int>>& A, const std::vector<int>& B ) const 
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

std::vector<complex_t> densityMatrix(const unsigned& dim, const std::vector<MinuitProxy>& pv )
{
  if( dim != 2 && dim != 3 )
  {
    std::vector<complex_t> rt(dim*dim);
    for( unsigned i = 0 ; i != dim; ++i ) rt[ dim*i +i ] = 1;
    return rt;
  }
  double px = pv[0];
  double py = pv[1];
  double pz = pv[2];
  if( dim == 2 ) return {1+pz    ,  px+1i*py,
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


std::vector<TransitionMatrix<void>> PolarisedSum::matrixElements() const
{
  return m_matrixElements;  
}


void   PolarisedSum::prepare()
{
  auto resetFlags  = [](auto& t){ t.workToDo = false; t.resetExternals() ; };
  auto flagUpdate  = [this](auto& t){ t.workToDo = this->m_nCalls == 0 || t.hasExternalsChanged(); };
  auto updateData  = [this](auto& t) mutable { if( t.workToDo && this->m_events != nullptr ) this->m_cache.update(this->m_events->store(), t) ; };
  auto updateInteg = [this](auto& t) mutable { if( t.workToDo ) this->m_integrator.updateCache(t) ; };

  transferParameters();
  for_each_sequence(m_matrixElements.begin(), m_matrixElements.end(), flagUpdate, updateData, updateInteg); 
  if( m_integrator.isReady() ) updateNorms();
  std::for_each( m_matrixElements.begin(), m_matrixElements.end(), resetFlags );
  if constexpr( detail::debug_type<PolarisedSum>::value )
  if( m_nCalls % 10000 == 0 ) debug_norm();
  m_pdfCache.update(m_cache, m_probExpression); 
  DEBUG( "m_pdfCache[0] = " << utils::at(m_pdfCache[0],0) << " w/o caching = " << getValNoCache(m_events->at(0)) << " w = " << m_weight << " N = " << m_norm );
  m_nCalls++; 
}

float_v PolarisedSum::operator()( const float_v*, const unsigned index ) const 
{ 
  return ( m_weight / m_norm ) * m_pdfCache[index];
}

#if ENABLE_AVX
double PolarisedSum::operator()( const double*, const unsigned index ) const 
{ 
  return operator()((const float_v*)nullptr, index / utils::size<float_v>::value ).at( index % utils::size<float_v>::value );
}
#endif


void PolarisedSum::debug_norm()
{
  if( !m_integrator.isReady() ) return; 
  double norm_slow = 0;
  for( auto& evt : *m_integrator.events<EventList_type>() ) 
    norm_slow += evt.weight() * getValNoCache(evt) / evt.genPdf();
  norm_slow /= m_integrator.norm();
  INFO("Norm: " << std::setprecision(10) << "bilinears=" << m_norm << "; Slow="   << norm_slow << "; d = "     << m_norm - norm_slow);
}

void   PolarisedSum::setEvents( EventList_type& events )
{ 
  reset();
  if( m_events != nullptr && m_ownEvents ) delete m_events; 
  m_events = &events;
  m_cache    . allocate( m_events->size(), m_matrixElements, m_dim.first * m_dim.second );
  m_pdfCache . allocate( m_events->size(), m_probExpression);
}

void   PolarisedSum::setMC( EventList_type& events )
{
  m_nCalls = 0;
  m_integrator = Integrator(&events, m_matrixElements, m_dim.first * m_dim.second );
  m_integIndex.clear();
  for( auto& i : m_matrixElements )
    m_integIndex.push_back( m_integrator.getCacheIndex(i) );
}

size_t PolarisedSum::size() const 
{ 
  return m_dim.first * m_dim.second * m_matrixElements.size(); 
}

void   PolarisedSum::reset( const bool& flag ){ m_nCalls = 0 ; }

Tensor PolarisedSum::transitionMatrix() const 
{ 
  auto size = m_dim.first * m_dim.second;
  std::vector<Expression> expressions(size, 0);
  unsigned totalSize = 0 ; 
  for( const auto& me : m_matrixElements ){
    auto coupling   = me.coupling.to_expression();
    INFO( me.decayDescriptor() << " " << coupling );
    auto cacheIndex = totalSize; 
    for( size_t i = 0 ; i < size ; ++i ){
      expressions[i] = expressions[i] + coupling * Parameter( "x1["+std::to_string(cacheIndex+i)+"]",0,true); 
    }
    totalSize += size; 
  }
  Tensor T_matrix(expressions, {m_dim.first, m_dim.second});
  T_matrix.st();
  return T_matrix; 
}

real_t PolarisedSum::operator()(const Event& evt) const 
{   
  return (m_weight/m_norm) * utils::at( m_pdfCache[ evt.index() / utils::size<float_v>::value ], evt.index() % utils::size<float_v>::value );
}

double PolarisedSum::norm() const 
{
  return m_norm;
}

complex_t PolarisedSum::norm(const size_t& i, const size_t& j, Integrator* integ)
{
  auto   ai = m_integIndex[i];
  auto   aj = m_integIndex[j];
  complex_t total = 0;
  auto    s1 = m_dim.first;
  auto    s2 = m_dim.second;
  for(size_t x = 0 ; x < m_norms.size(); ++x){
    auto f         = x % s2;
    auto psiIndex  = (x-f) / s2;
    auto m2        = psiIndex % s1;
    auto m1        = (psiIndex-m2)/s1;
    total          += m_rho[psiIndex] * m_norms[x].get(i, j, integ, ai+m1*s2+f, aj+m2*s2+f);
  }
  return total; 
}

void PolarisedSum::updateNorms()
{
  if(std::any_of(m_matrixElements.begin(),m_matrixElements.end(), [](auto& me){ return me.workToDo; } )){
  for( unsigned i = 0 ; i < m_matrixElements.size(); ++i ){
    for( unsigned j = i; j < m_matrixElements.size(); ++j ){
      if( m_matrixElements[i].workToDo || m_matrixElements[j].workToDo ) norm(i, j, &m_integrator);
    }
  }
  m_integrator.flush();
  }
  complex_t z = 0;
  for(size_t i = 0; i < m_matrixElements.size(); ++i){
    for(size_t j = 0; j < m_matrixElements.size(); ++j){
       z += m_matrixElements[i].coupling()*std::conj(m_matrixElements[j].coupling()) * ( i > j ? std::conj(norm(j,i)) : norm(i,j) );
    }
  }
  m_norm = std::real(z); 
}

void PolarisedSum::debug(const Event& evt)
{
  auto tsize = m_dim.first * m_dim.second;   
  std::vector<complex_v> this_cache; 
  for(unsigned j = 0; j != m_matrixElements.size(); ++j)
  {
    for(unsigned i = 0 ; i != tsize; ++i ) this_cache.emplace_back( m_cache(evt.index() / utils::size<float_v>::value, j*tsize + i) );
    INFO( m_matrixElements[j].decayDescriptor() << " " << vectorToString(this_cache, " ") );
    if( m_debug ) m_matrixElements[j].debug( evt ); 
  }   
  if( m_debug ) m_probExpression.debug(this_cache.data() );
  INFO("P(x) = " << getValNoCache(evt) << " " << operator()((const float_v*)nullptr, evt.index() / utils::size<float_v>::value ) );
  INFO("Prod = [" << vectorToString(m_pVector , ", ") <<"]");
}

void PolarisedSum::generateSourceCode(const std::string& fname, const double& normalisation, bool add_mt)
{
  INFO("Generating sourceCode -> " << fname );
  std::ofstream stream( fname );
  size_t size = m_dim.first * m_dim.second; 
  CompilerWrapper().preamble( stream );
  Expression event = Parameter("x0",0,true); 
  std::vector<Expression> expressions(size);
  for( auto& p : m_matrixElements ){
    auto expr = CompiledExpression<std::vector<complex_t>(const real_t*, const real_t*)>(
          p.expression(), 
          p.decayDescriptor(),
          m_eventType.getEventFormat(), DebugSymbols() ,m_mps ) ;
    expr.prepare();
    expr.to_stream( stream );
    expr.compileWithParameters( stream );
    Array z( make_cse( Function( programatic_name( p.name()) + "_wParams", {event} ) ), size );
    INFO( p.decayDescriptor() << " coupling = " << p.coupling() );
    for( unsigned int j = 0 ; j < size; ++j ) expressions[j] = expressions[j] + p.coupling() * z[j];
  }
  Tensor T_matrix( expressions, {m_dim.first, m_dim.second} );
  T_matrix.st();
  auto amp        = probExpression(T_matrix, convertProxies(m_pVector, [](auto& proxy) -> Expression{ return double(proxy);} )); 
  auto amp_extPol = probExpression(T_matrix, {Parameter("x2",0,true), Parameter("x3",0,true), Parameter("x4",0,true)}); 
  stream << CompiledExpression<double(
                               const double*, 
                               const int&)>( amp / normalisation, "FCN", m_mps, disableBatch() ) << std::endl ;

  stream << CompiledExpression<double(
                               const double*, 
                               const int&, 
                               const double&, 
                               const double&, 
                               const double&)>(amp_extPol / normalisation, "FCN_extPol", m_mps, disableBatch() ) << std::endl;
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
  ADD_DEBUG_TENSOR(T_conj, db);
  ADD_DEBUG_TENSOR(rho, db);
  ADD_DEBUG_TENSOR(TT , db);
  Expression rt = rho(a,b) * TT(b,a);
  return Real(rt);  
}

std::vector<FitFraction> PolarisedSum::fitFractions(const LinearErrorPropagator& prop)
{
  bool recomputeIntegrals    = NamedParameter<bool>("PolarisedSum::RecomputeIntegrals", false );
  bool interferenceFractions = NamedParameter<bool>("PolarisedSum::InterferenceFractions", false );
  std::vector<FitFraction> outputFractions; 
  for(const auto& rule : m_rules.rules()) 
  {
    FitFractionCalculator<PolarisedSum> pCalc(this, findIndices(m_matrixElements, rule.first), recomputeIntegrals);
    for(const auto& process : rule.second) 
    {
      if(process.head() == m_eventType.mother() && process.prefix() != m_prefix) continue;
      auto numeratorIndices   = processIndex(m_matrixElements, process.name());
      if(numeratorIndices.size() == 0 || numeratorIndices == pCalc.normSet ) continue; 
      pCalc.emplace_back(process.name(), numeratorIndices);
    }
    if( pCalc.calculators.size() == 0 ) continue;  
    auto fractions = pCalc(rule.first, prop);
    for( const auto& f : fractions ) outputFractions.emplace_back(f);
  }
  INFO("Fit fractions: ");
  for(const auto& p : outputFractions) INFO(p);  
  
  if( interferenceFractions )
  {  
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
    INFO("Interference fractions: ");
    for( auto& f : ifractions ) INFO( FitFraction(f) );
  }
  INFO("Returning: " << outputFractions.size() << " fractions");
  return outputFractions;
}

void PolarisedSum::transferParameters()
{ 
  m_probExpression.prepare();
  for(auto& me : m_matrixElements){
    me.coefficient = me.coupling();
    me.prepare();
  }
  for(auto& p : m_pVector) p.update();
  m_weight.update();
  m_rho = densityMatrix(m_dim.first, m_pVector);
}

real_t PolarisedSum::getValNoCache( const Event& evt ) const
{
  auto tsize = m_dim.first * m_dim.second;   
  std::vector<complex_v> cache( tsize * m_matrixElements.size() );
  for( unsigned i = 0 ; i != m_matrixElements.size(); ++i ){
    std::memmove( cache.data() + tsize * i , m_matrixElements[i](evt).data(), tsize * sizeof(complex_v) );
  }
  return utils::get<0>(m_probExpression( cache.data() ));
}

void   PolarisedSum::setWeight( MinuitProxy param ){ m_weight = param; } 
double PolarisedSum::getWeight() const { return m_weight ; }


std::function<real_t(const Event&)> PolarisedSum::evaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  Store<complex_v, Alignment::AoS> store(events->size(), m_matrixElements, m_dim.first * m_dim.second);
  for( auto& me : m_matrixElements ) store.update(events->store(), me );

  std::vector<double> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    utils::store(values.data() + utils::size<float_v>::value * block, (m_weight/m_norm) * m_probExpression(&store(block,0)) ); 
  }
  for( unsigned int i = 0 ; i != 10; ++i )
    DEBUG(values[i] << " " << getValNoCache( events->at(i) )  * ( m_weight / m_norm ) ); 
  return arrayToFunctor<double, Event>(values);
}

KeyedFunctors<double, Event> PolarisedSum::componentEvaluator(const EventList_type* events) const 
{
  auto& cache = m_integrator.cache();
  KeyedFunctors<double, Event> rt;
  for( unsigned i = 0 ; i != m_matrixElements.size(); ++i )
  {
    for( unsigned j = i ; j != m_matrixElements.size(); ++j ){
      auto mi = m_matrixElements[i]; 
      auto mj = m_matrixElements[j]; 
      auto ci = this->m_matrixElements[i].coefficient;
      auto cj = this->m_matrixElements[j].coefficient;
      double s = (i==j) ? 1 : 2 ;
      auto name = programatic_name(mi.decayTree.decayDescriptor()) + "_" + programatic_name( mj.decayTree.decayDescriptor() );
      INFO("Adding evaluator for: " << name  );
      auto functor = [ci,cj,i,j,s, &cache, this](const Event& event){  
      auto [s1,s2] = this->m_dim;
      auto R = s1 * s2; 
      complex_t total = 0;
      for( unsigned x = 0; x != this->m_norms.size(); ++x )
      {
        auto f         = x % s2;
        auto psiIndex  = (x-f) / s2;
        auto m2        = psiIndex % s1;
        auto m1        = (psiIndex-m2)/s1;
        total          += this->m_rho[psiIndex] * ci * cache.get<complex_t>(event.index(),R * i + m1 * s2 + f) 
                                     * std::conj( cj * cache.get<complex_t>(event.index(),R * j + m2 * s2 + f) );
      }
      return s * std::real(total);
      };
      rt.add(functor, name, "");
    }
  }
  return rt; 
}
