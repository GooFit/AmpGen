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
  : m_rules    (mps)
  , m_evtType  (type)
  , m_printFreq(NamedParameter<size_t>(     "CoherentSum::PrintFrequency", 100)  )
  , m_dbThis   (NamedParameter<bool>(       "CoherentSum::Debug"         , false))
  , m_verbosity(NamedParameter<bool>(       "CoherentSum::Verbosity"     , 0)    )
  , m_objCache (NamedParameter<std::string>("CoherentSum::ObjectCache"   ,"")    )
  , m_prefix   (prefix)
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
    m_matrixElements[i] = TransitionMatrix<CoherentSum::complex_v>( amplitudes[i].first, amplitudes[i].second, mps, this->m_evtType.getEventFormat(), this->m_dbThis);
    CompilerWrapper().compile( m_matrixElements[i].amp, this->m_objCache); } );
  }
  m_isConstant = false;
}


void CoherentSum::prepare()
{
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  transferParameters(); 
  std::vector<size_t> changedPdfIndices;
  ProfileClock clockEval; 
  if( m_prepareCalls == 0 && m_events != nullptr ){
    m_events->reserveCache(m_matrixElements.size());
    for( auto& me : m_matrixElements ) me.addressData = m_events->registerExpression( me.amp );
  }
  if( m_prepareCalls == 0 ) m_integrator.allocate( m_matrixElements );
  for ( size_t i = 0; i < m_matrixElements.size(); ++i ) {
    m_matrixElements[i].amp.prepare();
    if ( m_prepareCalls != 0 && !m_matrixElements[i].amp.hasExternalsChanged() ) continue;
    if ( m_events != nullptr ) m_events->updateCache( m_matrixElements[i].amp, m_matrixElements[i].addressData ); 
    m_integrator.prepareExpression(m_matrixElements[i].amp );
    changedPdfIndices.push_back(i);
    m_matrixElements[i].amp.resetExternals();
  }
  clockEval.stop();
  ProfileClock clockIntegral;
  if ( m_integrator.isReady())  updateNorms( changedPdfIndices );
  else if ( m_verbosity ) WARNING( "No simulated sample specified for " << this );
  m_norm = norm();
  if ( m_prepareCalls == 0 ){
    INFO( "Norm: " << m_norm );
    for(unsigned i = 0 ; i != m_matrixElements.size() ; ++i ){
      for(unsigned j = 0 ; j != m_matrixElements.size() ; ++j ){
        if( std::isnan( std::real(m_normalisations(i,j) )) || std::isnan( std::imag(m_normalisations(i,j))) ) 
          ERROR("Norm: " << m_matrixElements[i].name() << " " << m_matrixElements[j].name() << " is ill-posed!");
      }
    }
    // INFO( m_normalisations.get(0,0) << " "  
   //    << m_normalisations.get(1,0) << " "  
   //    << m_normalisations.get(0,1) << " "  
   //    << m_normalisations.get(2,2) ); 
  }
  if ( m_verbosity && changedPdfIndices.size() !=0 ) {
    clockIntegral.stop();
    INFO( "Time Performance: "
        << "Eval = "       << clockEval     << " ms"
        << ", Integral = " << clockIntegral << " ms"
        << ", Total = "    << clockEval + clockIntegral << " ms; normalisation = "  << m_norm );
    m_lastPrint = m_prepareCalls;
  }
  m_prepareCalls++;
}

void CoherentSum::updateNorms( const std::vector<size_t>& changedPdfIndices )
{
  std::vector<size_t> cacheIndex;
  std::transform( m_matrixElements.begin(), m_matrixElements.end(), std::back_inserter(cacheIndex), 
    [this](auto& m){ return this->m_integrator.getCacheIndex( m.amp ) ; } );
  for ( auto& i : changedPdfIndices )
    for ( size_t j = 0; j < size(); ++j )
      m_integrator.queueIntegral( cacheIndex[i], cacheIndex[j] ,i, j, &m_normalisations );
  m_integrator.flush();
  m_normalisations.resetCalculateFlags();
}

void CoherentSum::debug( const Event& evt, const std::string& nameMustContain )
{
  prepare();
  for ( auto& me : m_matrixElements ) me.amp.resetExternals();
  if ( nameMustContain == "" )
    for ( auto& me : m_matrixElements ) {
      auto A = me(evt);
      INFO( std::setw(70) << me.decayTree.uniqueString() 
          << " A = [ "  << utils::get(A.real())             << " " << utils::get(A.imag())
          << " ] g = [ "<< me.coupling().real() << " " << me.coupling().imag() << " ] "
          << m_events->cache( evt.index(), me.addressData )
          << me.decayTree.CP() );

     // if( m_dbThis ) me.amp.debug( evt.address() );
    }
  //else
  //  for ( auto& me : m_matrixElements )
  //    if ( me.amp.name().find( nameMustContain ) != std::string::npos ) me.amp.debug( evt.address() );
  // if( evt.cacheSize() != 0 ) INFO( "Pdf = " << prob_unnormalised( evt ) );
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
    stream << p.amp << std::endl;
    p.amp.compileWithParameters( stream );
    if( includePythonBindings ) p.amp.compileDetails( stream );
  }
  Expression event = Parameter("x0",0,true);
  Expression pa    = Parameter("double(x1)",0,true);
  Expression amplitude;
  for( unsigned int i = 0 ; i < size(); ++i ){
    auto& p = m_matrixElements[i];
    Expression this_amplitude = p.coupling() * Function( programatic_name( p.amp.name() ) + "_wParams", {event} ); 
    amplitude = amplitude + ( p.decayTree.finalStateParity() == 1 ? 1 : pa ) * this_amplitude; 
  }
  stream << CompiledExpression<std::complex<double>(const double*, const int&)>( amplitude  , "AMP" ) << std::endl; 
  stream << CompiledExpression<double(const double*, const int&)>(fcn::norm(amplitude) / normalisation, "FCN" ) << std::endl; 
  if( includePythonBindings ){
    stream << CompiledExpression<unsigned int(void)>( m_matrixElements.size(), "matrix_elements_n" ) << std::endl;
    stream << CompiledExpression<double      (void)>( normalisation, "normalization") << std::endl;

    stream << "extern \"C\" const char* matrix_elements(int n) {\n";
    for ( size_t i = 0; i < m_matrixElements.size(); i++ ) {
      stream << "  if(n ==" << i << ") return \"" << m_matrixElements.at(i).amp.progName() << "\" ;\n";
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
      stream << programatic_name( p.amp.name() )<< "_wParams( E )";
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
  for ( auto& mE : m_matrixElements ) mE.addressData = 999;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integrator = Integrator_type();
  }
}

void CoherentSum::setEvents( EventList_type& list )
{
  if ( m_verbosity ) INFO( "Setting event list with:" << list.size() << " events for " << this );
  reset();
  m_events = &list;
}


void CoherentSum::setMC( EventList_type& sim )
{
  if ( m_verbosity ) INFO( "Setting norm. event list with:" << sim.size() << " events for " << this );
  reset();
  m_integrator = Integrator_type( &sim );
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
      auto val = norms.get(i, j) 
        * m_matrixElements[i].coefficient 
        * std::conj(m_matrixElements[j].coefficient);
      acc += val;
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
    unsigned int address = mE.addressData;
    std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << m_events->cache( evt.index(), address )
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
  for ( const auto& mE : m_matrixElements ) {
    value = value + mE.coefficient * m_events->cache( evt.index(), mE.addressData );
  }
  #if ENABLE_AVX2 
  return value.at(evt.index() % float_v::size);
  #else 
  return value;
  #endif
}

#if ENABLE_AVX2 
float_v CoherentSum::operator()( const float_v* /*evt*/, const unsigned block ) const 
{
  complex_v value( 0., 0. );
  for ( const auto& mE : m_matrixElements ) {
    value = value + mE.coefficient * m_events->cache()[ block * m_events->cacheSize() + mE.addressData ];
  }
  return (m_weight/m_norm ) * AVX2::norm( value ); 
}

#endif


std::function<real_t(const Event&)> CoherentSum::evaluator(const EventList_type* events) const 
{
  if( events != nullptr && events != m_integrator.events() )
    ERROR("Evaluator only working on the integration sample, fix me!"); 
  std::vector<unsigned> address_mapping( size() );
  for( const auto& me : m_matrixElements ) address_mapping[me.addressData] = m_integrator.getCacheIndex( me.amp );
  std::vector<double> values( events->size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int i = 0 ; i != events->size(); ++i )
  {
    complex_t amp = 0;
    for( unsigned j = 0 ; j != address_mapping.size(); ++j ) amp += m_matrixElements[j].coefficient * this->m_integrator.get(address_mapping[j], i);
    values[i] = m_weight * std::norm(amp) / m_norm;
  }
  return arrayToFunctor<double, Event>(values);
}

KeyedView<double, CoherentSum::EventList_type> CoherentSum::componentEvaluator(const EventList_type* events) const 
{
  if( events != nullptr && events != m_integrator.events() ) 
    ERROR("Evaluator only working on the integration sample, fix me!"); 
  
  KeyedView<double, EventList_type> rt(*events, m_matrixElements.size() ); 
  std::vector<unsigned> address_mapping(m_matrixElements.size());
  for( unsigned i = 0; i != m_matrixElements.size(); ++i ) address_mapping[i] = m_integrator.getCacheIndex( m_matrixElements[i].amp ); 

  for( unsigned i = 0 ; i != m_matrixElements.size(); ++i )
  {
    auto& me = m_matrixElements[i]; 
    rt.setKey(i, programatic_name( me.decayTree.decayDescriptor() ) );
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( unsigned evt = 0 ; evt != events->size(); ++evt )
    {
      complex_t total = 0;
      for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ){
          total          +=  this->m_integrator.get( address_mapping[i], evt ) * m_matrixElements[i].coefficient 
                * std::conj( this->m_integrator.get( address_mapping[j], evt ) * m_matrixElements[j].coefficient );
      } 
      rt(events->at(evt), i) = m_weight * std::real( total ) / m_norm;  
    }
  }
  return rt; 
}

