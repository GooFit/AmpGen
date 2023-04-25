#include "AmpGen/CoherentSum.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ratio>
#include <thread>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/ThreadPool.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/simd/utils.h"
#include "AmpGen/Array.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AmpGen;
CoherentSum::CoherentSum() = default; 

CoherentSum::CoherentSum( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix )
  :   m_evtType  (type)
      , m_printFreq(NamedParameter<size_t>(     "CoherentSum::PrintFrequency", 100)  )
      , m_dbThis   (NamedParameter<bool>(       "CoherentSum::Debug"         , false))
      , m_verbosity(NamedParameter<bool>(       "CoherentSum::Verbosity"     , 0)    )
  , m_objCache (NamedParameter<std::string>("CoherentSum::ObjectCache"   ,"")    )
  , m_prefix   (prefix)
      , m_mps(&mps) 
{
  auto rules = AmplitudeRules::create(mps);
  auto amplitudes      = rules->getMatchingRules( m_evtType, prefix);
  if( amplitudes.size() == 0 ){
    WARNING("The defined amplitudes don't seem to be able to be able to generate eventType: " << type);
  }
  for( auto& amp : amplitudes ) INFO( prefix + amp.first.decayDescriptor() );
  m_matrixElements.resize( amplitudes.size() );
  m_normalisations.resize( m_matrixElements.size(), m_matrixElements.size() ); 
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , std::thread::hardware_concurrency(), "Number of threads to use" );
  ThreadPool tp(nThreads);
  auto head_rules = rules->rulesForDecay(m_evtType.mother(), m_prefix);
  /*
  for( const auto& rule : head_rules )
  {
    auto indices = findIndices(m_matrixElements, rule.particle().decayDescriptor() );
    INFO( rule.particle() << " [" << vectorToString(indices, " ") << "]" ); 
  }
  */
  for(size_t i = 0; i < m_matrixElements.size(); ++i){
    tp.enqueue( [i, this, &mps, &amplitudes]() mutable {
        this->m_matrixElements[i] = 
          MatrixElement(amplitudes[i].first, amplitudes[i].second, mps, this->m_evtType.getEventFormat(), this->m_dbThis);  
        CompilerWrapper().compile( this->m_matrixElements[i], this->m_objCache); 
    } ); 
  }
}

void CoherentSum::prepare()
{
  transferParameters(); 
  ProfileClock clockEval; 
  for (auto& t : m_matrixElements ) {
    if( not t.isReady() )
      FATAL( t.decayDescriptor() << " not ready, fix me"); 
    t.prepare();
    
    if ( m_prepareCalls != 0 && !t.hasExternalsChanged() ) continue;
    if ( m_events != nullptr ) m_cache.update(t);
    m_integrator.updateCache(t);
    t.resetExternals();
    t.workToDo = true; 
  }
  clockEval.stop();
  ProfileClock clockIntegral;
  if (m_integrator.isReady()) updateNorms();
  else if ( m_verbosity ) WARNING( "No simulated sample specified for " << this );
  clockIntegral.stop();
  if ( m_verbosity && m_prepareCalls % 100 == 0  ) {
    INFO( "Time Performance: "
        << "Eval = "       << clockEval     << " ms"
        << ", Integral = " << clockIntegral << " ms"
        << ", Total = "    << clockEval + clockIntegral << " ms; normalisation = "  << m_norm );
    m_lastPrint = m_prepareCalls;
    
  }
  // for( int i = 0 ; i != m_cache.nFields(); ++i ) 
  //   INFO( m_matrixElements[i].name() << " " <<  std::setprecision(15) << m_cache(0, m_cache.find(m_matrixElements[i].name())[0] ));
  // INFO( (*this)( nullptr, 0 )  << " " << m_norm ); 
  for( auto& t : m_matrixElements ) t.workToDo = false; 
  m_prepareCalls++;
}

void CoherentSum::updateNorms()
{
  if(std::any_of(m_matrixElements.begin(),m_matrixElements.end(), [](auto& me){ return me.workToDo; } ))
  {
    for ( unsigned i = 0; i != m_matrixElements.size(); ++i )
    {
      for ( size_t j = i; j < size(); ++j ){
        if( m_matrixElements[i].workToDo || m_matrixElements[j].workToDo ) 
          m_normalisations.get(i, j, &m_integrator, i, j);
      }
    }
  }
  m_integrator.flush();
  m_normalisations.resetCalculateFlags();
  m_norm = norm();
}

void CoherentSum::debug( const Event& evt, const std::string& nameMustContain )
{
  prepare();
  INFO("Weight = " << evt.weight() << " genPDF = " << evt.genPdf() );

  for ( auto& me : m_matrixElements ) {
    auto A = me(evt);
    INFO( std::setw(70) << me.decayTree.uniqueString() 
        << " A = [ "  << A[0].real()             << " " << A[0].imag()
        << " ] g = [ "<< me.coupling().real() << " " << me.coupling().imag() << " ] "
        << m_cache( evt.index(), std::distance(&m_matrixElements[0], &me ) )
        << me.decayTree.CP() );
  }
  if( m_dbThis ) for ( auto& me : m_matrixElements ) me.debug( evt) ;
  INFO( "A(x) = " << getVal(evt) << " without cache: " << getValNoCache(evt) );
}

std::vector<FitFraction> CoherentSum::fitFractions(const LinearErrorPropagator& linProp)
{
  prepare();
  bool recomputeIntegrals    = NamedParameter<bool>("CoherentSum::RecomputeIntegrals", false );
  std::vector<FitFraction> outputFractions;
  auto rules = AmplitudeRules::get();
  for(auto& rule : rules->rules() ) 
  {
    FitFractionCalculator<CoherentSum> pCalc(this, findIndices(m_matrixElements, rule.first), recomputeIntegrals);
    for(auto& process : rule.second) 
    {
      if(process.head() == m_evtType.mother() && process.prefix() != m_prefix) continue;
      auto numeratorIndices   = processIndex(m_matrixElements,process.name());
      if(numeratorIndices.size() == 0 || numeratorIndices == pCalc.normSet ) continue; 
      pCalc.emplace_back(process.name(), numeratorIndices);
    }
    if( pCalc.calculators.size() == 0 ) continue;  
    auto fractions = pCalc(rule.first, linProp);
    std::transform( fractions.begin(), fractions.end(), std::back_inserter(outputFractions),[](auto& p){ return p;} );
  };
  auto ffForHead = rules->rulesForDecay(m_evtType.mother(), m_prefix);
  FitFractionCalculator<CoherentSum> iCalc(this, findIndices(m_matrixElements, m_evtType.mother()), recomputeIntegrals);
  for ( size_t i = 0; i < ffForHead.size(); ++i ) 
  {
    for ( size_t j = i + 1; j < ffForHead.size(); ++j ) 
    {
      iCalc.emplace_back(ffForHead[i].name() + "x" + ffForHead[j].name(), 
          processIndex(m_matrixElements, ffForHead[i].name()), 
          processIndex(m_matrixElements, ffForHead[j].name()) );
    }
  }
  std::vector<FitFraction> interferenceFractions = iCalc(m_evtType.mother()+"_interference",linProp);
  for(auto& f : interferenceFractions) outputFractions.push_back(f); 
  for(auto& p : outputFractions) INFO(p);
  return outputFractions;
}

void CoherentSum::generateSourceCode(const std::string& fname, const double& normalisation, bool add_mt)
{
  std::vector<CompiledExpressionBase*> functions; 
  transferParameters();
  Tensor rt( std::vector<unsigned>{unsigned(m_matrixElements.size())} ); 
  Expression event = Parameter("x0",0,true);
  Expression pa    = Parameter("double(x1)",0,true);
  Expression amplitude;
  for( unsigned i = 0 ; i != size(); ++i ) rt[i] = m_matrixElements[i].expression(); 
  functions.push_back( new CompiledExpression<std::vector<complex_t>( const real_t*, const real_t* )>(TensorExpression(rt), "all_amplitudes", m_evtType.getEventFormat(), includeParameters(), m_mps, disableBatch() ) );  
  auto ampTensor = Array( make_cse(Function("all_amplitudes_wParams", {event} )), -1 ); 
  for( unsigned int i = 0 ; i < size(); ++i ){
    auto& p = m_matrixElements[i];
    Expression this_amplitude = p.coupling() * ampTensor[i];
    amplitude += ( p.decayTree.finalStateParity() == 1 ? 1 : pa ) * this_amplitude; 
    functions.push_back( new CompiledExpression<complex_t( const real_t*)>( ampTensor[i], p.decayDescriptor() + "_wParams", m_evtType.getEventFormat(), disableBatch() ) );  
  }
  functions.push_back( new CompiledExpression<complex_t(const real_t*, const int&)>( amplitude  , "AMP", disableBatch() ) );
  functions.push_back( new CompiledExpression<real_t(const real_t*, const int&)>(fcn::norm(amplitude) / normalisation, "FCN", disableBatch() ) ); 

  CompilerWrapper().compile( functions, fname );

  for( auto& function : functions ) delete function;  

  INFO("Generating source-code for PDF: " << fname << " include MT symbols? " << add_mt << " normalisation = " << normalisation );
}

complex_t CoherentSum::getValNoCache( const Event& evt ) const
{
  auto v = complex_v(std::accumulate( m_matrixElements.begin(), 
        m_matrixElements.end(), 
        complex_v(0,0), 
        [&evt]( const auto& a, const auto& b ){ return a + complex_v(b.coefficient) * b(evt)[0];} ));
  return complex_t( utils::get<0>(v.real()), utils::get<0>(v.imag()) ); 
}

void CoherentSum::reset( bool resetEvents )
{
  m_prepareCalls                                     = 0;
  m_lastPrint                                        = 0;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integrator = Integrator();
  }
}

void CoherentSum::setEvents( const EventList_type& list )
{
  DEBUG( "Setting event list with:" << list.size() << " events for " << this );
  reset();
  for( auto& me : m_matrixElements ){ DEBUG("Registering: " << me.name() ) ; }
  if( m_ownEvents && m_events != nullptr ) delete m_events; 
  m_events = &list;
  m_cache.allocate(m_events, m_matrixElements); 
}


void CoherentSum::setMC( const EventList_type& sim )
{
  if ( m_verbosity ) INFO( "Setting norm. event list with:" << sim.size() << " events for " << this );
  reset();
  m_integrator = Integrator( &sim, m_matrixElements );
}

real_t CoherentSum::norm() const
{
  return norm(m_normalisations);
}

real_t CoherentSum::norm(const Bilinears& norms) const
{
  complex_t acc(0, 0);
  for ( size_t i = 0; i < size(); ++i ) {
    for ( size_t j = 0; j < size(); ++j ) {
      // INFO( i << " " << j << " " << m_matrixElements[i].coefficient * std::conj(m_matrixElements[j].coefficient) << " " <<  ( i > j ? std::conj(norm(j,i)) : norm(i,j) ) );
      acc += m_matrixElements[i].coefficient * std::conj(m_matrixElements[j].coefficient)* ( i > j ? std::conj(norm(j,i)) : norm(i,j) );
    }
  }
  return acc.real();
}

complex_t CoherentSum::norm(const size_t& x, const size_t& y) const
{
  return m_normalisations.get(x, y);
}

void CoherentSum::transferParameters()
{
  for ( auto& mE : m_matrixElements ) mE.coefficient = mE.coupling();
  m_weight.update();
}

void CoherentSum::printVal(const Event& evt)
{
  /*
     for ( auto& mE : m_matrixElements ) {
     unsigned int address = std::distance( &mE , &m_matrixElements[0] );
     std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << m_cache( evt.index() / utils::size<real_v>::value, address )
     << " address = " << address << " " << mE( evt ) << std::endl;
     if( mE.coupling.size() != 1 ){
     std::cout << "CouplingConstants: " << std::endl;
     mE.coupling.print();
     std::cout << "================================" << std::endl;
     }
     }
     */
}

complex_t CoherentSum::getVal( const Event& evt ) const
{
  complex_v value( 0., 0. );
  for (unsigned int i = 0 ; i != m_matrixElements.size(); ++i ) {
    value += complex_v( m_matrixElements[i].coefficient ) * m_cache(evt.index() / utils::size<real_v>::value, i );
  }
#if ENABLE_AVX
  return utils::at(value, evt.index() % utils::size<real_v>::value);
#else 
  return value;
#endif
}

real_v CoherentSum::operator()( const real_v* /*evt*/, const unsigned block ) const 
{
  complex_v value( 0., 0. );
  for ( const auto& mE : m_matrixElements ) 
  {
    unsigned address = &mE - &m_matrixElements[0];
    value += complex_v(mE.coefficient) * m_cache(block, address); 
  }
  return (m_weight/m_norm ) * utils::norm(value);
}

#if ENABLE_AVX

//double CoherentSum::operator()( const double* /*evt*/, const unsigned block ) const 
//{
//  return operator()((const real_v*)nullptr, block / utils::size<real_v>::value ).at( block % utils::size<real_v>::value );
//}
#endif

std::function<real_t(const Event&)> CoherentSum::evaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  FunctionCache<EventList_type, complex_v, Alignment::AoS> store(events, m_matrixElements);
  for( auto& me : m_matrixElements ) store.update(me);

  std::vector<double> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    complex_v amp(0.,0.);
    for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ) 
      amp = amp + complex_v(m_matrixElements[j].coefficient) * store(block, j);
    utils::store( values.data() + block * utils::size<real_v>::value,  (m_weight/m_norm) * utils::norm(amp)  );
  }
  return arrayToFunctor<real_t, typename EventList_type::value_type>(values);
}

std::function<complex_t(const Event&)> CoherentSum::amplitudeEvaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  FunctionCache<EventList_type, complex_v, Alignment::AoS> store(events, m_matrixElements);
  for( auto& me : m_matrixElements ) store.update( me );
  std::vector<complex_t> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    complex_v amp(0.,0.);
    for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ) 
      amp = amp + complex_v(m_matrixElements[j].coefficient) * store(block, j);
    for( unsigned k = 0; k != utils::size<complex_v>::value; ++k )
    {
      values[ block * utils::size<complex_v>::value + k] = utils::at( amp, k ); 
    }
  }
  return arrayToFunctor<complex_t, typename EventList_type::value_type>(values);
}



KeyedFunctors<double(Event)> CoherentSum::componentEvaluator(const EventList_type* ievents) const 
{
  using store_t = FunctionCache<EventList_type, complex_v, Alignment::SoA>; 
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  KeyedFunctors<double(Event)> rt; 
  std::shared_ptr<const store_t> cache;
  if( events != m_integrator.events<EventList_type>() )
  {
    cache = std::make_shared<const store_t>(events, m_matrixElements);
    for( auto& me : m_matrixElements ) const_cast<store_t*>(cache.get())->update(me);
  }
  else cache = std::shared_ptr<const store_t>( & m_integrator.cache(), [](const store_t* t){} ); 
  /// this little slice of weirdness allows either a new cache to be instantiated, or one to just get a pointer to the one used for the integration. 

  for( unsigned i = 0 ; i != m_matrixElements.size(); ++i )
  {
    for( unsigned j = i ; j != m_matrixElements.size(); ++j )
    {
      auto mi = m_matrixElements[i]; 
      auto mj = m_matrixElements[j]; 
      auto ci = this->m_matrixElements[i].coefficient;
      auto cj = this->m_matrixElements[j].coefficient;
      double s = (i==j) ? 1 : 2 ;
      auto name = programatic_name(mi.decayTree.decayDescriptor()) + "_" + programatic_name( mj.decayTree.decayDescriptor() );
      auto functor = [ci,cj,i,j,s, cache](const Event& event){ return s * std::real( ci * cache->get<complex_t>( event.index(), i ) *  std::conj( cj * cache->get<complex_t>( event.index(), j ) ) ) ;};
      rt.add(functor, name, "");
    }
  }
  INFO("Returning ... "); 
  return rt; 
}

CoherentSum::~CoherentSum()
{
  if( m_ownEvents && m_events !=nullptr ) delete m_events; 
}
