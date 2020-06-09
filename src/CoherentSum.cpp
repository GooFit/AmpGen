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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AmpGen;
CoherentSum::CoherentSum() = default; 

CoherentSum::CoherentSum( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix )
  :   m_rules    (mps)
    , m_evtType  (type)
    , m_printFreq(NamedParameter<size_t>(     "CoherentSum::PrintFrequency", 100)  )
    , m_dbThis   (NamedParameter<bool>(       "CoherentSum::Debug"         , false))
    , m_verbosity(NamedParameter<bool>(       "CoherentSum::Verbosity"     , 0)    )
    , m_objCache (NamedParameter<std::string>("CoherentSum::ObjectCache"   ,"")    )
    , m_prefix   (prefix)
    , m_mps(&mps) 
{
  auto amplitudes      = m_rules.getMatchingRules( m_evtType, prefix);
  if( amplitudes.size() == 0 ){
    WARNING("The defined amplitudes don't seem to be able to be able to generate eventType: " << type);
  }
  for( auto& amp : amplitudes ) INFO( amp.first.decayDescriptor() );
  m_matrixElements.resize( amplitudes.size() );
  m_normalisations.resize( m_matrixElements.size(), m_matrixElements.size() ); 
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , std::thread::hardware_concurrency(), "Number of threads to use" );
  ThreadPool tp(nThreads);
  for(size_t i = 0; i < m_matrixElements.size(); ++i){
    tp.enqueue( [i,this,&mps,&amplitudes]{
       auto& [p, c] = amplitudes[i];
        DebugSymbols db;
        m_matrixElements[i] = 
          TransitionMatrix<complex_v>(p, c, 
          CompiledExpression<complex_v(const real_t*, const float_v*)>( p.getExpression(m_dbThis ? &db : nullptr), p.decayDescriptor(),
            this->m_evtType.getEventFormat(), db, &mps ) );
        CompilerWrapper().compile( m_matrixElements[i], this->m_objCache); 
      } ); 
  }
}

void CoherentSum::prepare()
{
  transferParameters(); 
  ProfileClock clockEval; 
  for (auto& t : m_matrixElements ) {
    t.prepare();
    if ( m_prepareCalls != 0 && !t.hasExternalsChanged() ) continue;
    if ( m_events != nullptr ) m_cache.update(m_events->store(), t );
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
  for( auto& t : m_matrixElements ) t.workToDo = false; 
  m_prepareCalls++;
}

void CoherentSum::updateNorms()
{
  if(std::any_of(m_matrixElements.begin(),m_matrixElements.end(), [](auto& me){ return me.workToDo; } ))
  {
    for ( unsigned i = 0; i != m_matrixElements.size(); ++i )
      for ( size_t j = i; j < size(); ++j ){
        if( m_matrixElements[i].workToDo || m_matrixElements[j].workToDo ) 
          m_normalisations.get(i, j, &m_integrator, i, j);
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
        << " A = [ "  << utils::get<0>(A.real())             << " " << utils::get<0>(A.imag())
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
  for(auto& rule : m_rules.rules()) 
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
  auto ffForHead = m_rules.rulesForDecay(m_evtType.mother(), m_prefix);
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
  std::ofstream stream( fname );
  transferParameters();
  stream << std::setprecision( 10 );
  stream << "#include <complex>\n";
  stream << "#include <vector>\n";
  stream << "#include <math.h>\n";
  if ( add_mt ) stream << "#include <thread>\n";
  bool includePythonBindings = NamedParameter<bool>("CoherentSum::IncludePythonBindings",false);

  for ( auto& p : m_matrixElements ){
    auto expr = CompiledExpression<complex_t(const real_t*, const real_t*)>(
          p.expression(), 
          p.decayDescriptor(),
          m_evtType.getEventFormat(), DebugSymbols() , m_mps );
    expr.prepare();
    stream << expr << std::endl;
    expr.compileWithParameters( stream );
    if( includePythonBindings ) p.compileDetails( stream );
  }
  Expression event = Parameter("x0",0,true);
  Expression pa    = Parameter("double(x1)",0,true);
  Expression amplitude;
  for( unsigned int i = 0 ; i < size(); ++i ){
    auto& p = m_matrixElements[i];
    Expression this_amplitude = p.coupling() * Function( programatic_name( p.name() ) + "_wParams", {event} ); 
    amplitude = amplitude + ( p.decayTree.finalStateParity() == 1 ? 1 : pa ) * this_amplitude; 
  }
  stream << CompiledExpression<std::complex<double>(const double*, const int&)>( amplitude  , "AMP" ) << std::endl; 
  stream << CompiledExpression<double(const double*, const int&)>(fcn::norm(amplitude) / normalisation, "FCN" ) << std::endl; 
  if( includePythonBindings ){
    stream << CompiledExpression<unsigned int(void)>( m_matrixElements.size(), "matrix_elements_n" ) << std::endl;
    stream << CompiledExpression<double      (void)>( normalisation, "normalization") << std::endl;

    stream << "extern \"C\" const char* matrix_elements(int n) {\n";
    for ( size_t i = 0; i < m_matrixElements.size(); i++ ) {
      stream << "  if(n ==" << i << ") return \"" << m_matrixElements.at(i).progName() << "\" ;\n";
    }
    stream << "  return 0;\n}\n";
    stream << "extern \"C\" void FCN_all(double* out, double* events, unsigned int size, int parity, double* amps){\n";
    stream << "  double local_amps [] = {\n";
    for ( unsigned int i = 0; i < size(); ++i ) {
      auto& p = m_matrixElements[i];
      stream << "    " << p.coupling().real() << ", " << p.coupling().imag() << ( i + 1 == size() ? "" : "," ) << "\n";
    }
    stream << "  };\n";
    stream << "  if(amps == nullptr)\n";
    stream << "    amps = local_amps;\n\n";
    stream << "  for(unsigned int i=0; i<size; i++) {\n";
    stream << "    double* E = events + i*" << ( *this->m_events )[0].size() << ";\n";

    stream << "  std::complex<double> amplitude = \n";
    for ( unsigned int i = 0; i < size(); ++i ) {
      stream << "    ";
      auto& p    = m_matrixElements[i];
      int parity = p.decayTree.finalStateParity();
      if ( parity == -1 ) stream << "double(parity) * ";
      stream << "std::complex<double>(amps[" << i * 2 << "],amps[" << i * 2 + 1 << "]) * ";
      stream << programatic_name( p.name() )<< "_wParams( E )";
      stream << ( i == size() - 1 ? ";" : " +" ) << "\n";
    }
    stream << "  out[i] =  std::norm(amplitude) / " << normalisation << ";\n  }\n}\n";

    stream << "extern \"C\" double coefficients( int n, int which, int parity){\n";
    for ( size_t i = 0; i < size(); i++ ) {
      auto& p    = m_matrixElements[i];
      int parity = p.decayTree.finalStateParity();
      stream << "  if(n == " << i << ") return ";
      if ( parity == -1 ) stream << "double(parity) * ";
      stream << "(which==0 ? " << p.coupling().real() << " : " << p.coupling().imag() << ");\n";
    }
    stream << "  return 0;\n}\n";
  }
  INFO("Generating source-code for PDF: " << fname << " include MT symbols? " << add_mt << " normalisation = " << normalisation );
  stream.close();
}

complex_t CoherentSum::getValNoCache( const Event& evt ) const
{
  return utils::get<0>( complex_v(std::accumulate( m_matrixElements.begin(), 
          m_matrixElements.end(), 
          complex_v(0,0), 
          [&evt]( const auto& a, const auto& b ){ return a + b.coefficient * b(evt);} )) );
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
  m_cache . allocate( m_events->size(), m_matrixElements );   
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
      acc += m_matrixElements[i].coefficient 
        * std::conj(m_matrixElements[j].coefficient)
        * ( i > j ? std::conj(norm(j,i)) : norm(i,j) );
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
  for ( auto& mE : m_matrixElements ) {
    unsigned int address = std::distance( &mE , &m_matrixElements[0] );
    std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << m_cache( evt.index() / utils::size<float_v>::value, address )
      << " address = " << address << " " << mE( evt ) << std::endl;
    if( mE.coupling.size() != 1 ){
      std::cout << "CouplingConstants: " << std::endl;
      mE.coupling.print();
      std::cout << "================================" << std::endl;
    }
  }
}

complex_t CoherentSum::getVal( const Event& evt ) const
{
  complex_v value( 0., 0. );
  for (unsigned int i = 0 ; i != m_matrixElements.size(); ++i ) {
    value = value + m_matrixElements[i].coefficient * m_cache(evt.index() / utils::size<float_v>::value, i );
  }
#if ENABLE_AVX
  return value.at(evt.index() % utils::size<float_v>::value );
#else 
  return value;
#endif
}

float_v CoherentSum::operator()( const float_v* /*evt*/, const unsigned block ) const 
{
  complex_v value( 0., 0. );
  for ( const auto& mE : m_matrixElements ) 
  {
    unsigned address = &mE - &m_matrixElements[0];
    value = value + mE.coefficient * m_cache(block, address); 
  }
  return (m_weight/m_norm ) * utils::norm(value); 
}

#if ENABLE_AVX
double CoherentSum::operator()( const double* /*evt*/, const unsigned block ) const 
{
  return operator()((const float_v*)nullptr, block / utils::size<float_v>::value ).at( block % utils::size<float_v>::value );
}
#endif

std::function<real_t(const Event&)> CoherentSum::evaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  Store<complex_v, Alignment::AoS> store( events->size(), m_matrixElements);
  for( auto& me : m_matrixElements ) store.update(events->store(), me );
  
  std::vector<double> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    complex_v amp(0.,0.);
    for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ) 
      amp = amp + m_matrixElements[j].coefficient * store(block, j);
    utils::store( values.data() + block * utils::size<float_v>::value,  (m_weight/m_norm) * utils::norm(amp)  );
  }
  return arrayToFunctor<double, typename EventList_type::value_type>(values);
}

KeyedFunctors<double, Event> CoherentSum::componentEvaluator(const EventList_type* ievents) const 
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
      auto functor = [ci,cj,i,j,s, &cache](const Event& event){ return s * std::real( ci * cache.get<complex_t>( event.index(), i ) *  std::conj( cj * cache.get<complex_t>( event.index(), j ) ) ) ;};
      rt.add(functor, name, "");
    }
  }
  INFO(" Returning: " << rt.keys.size() << " functors" );
  return rt; 
}
