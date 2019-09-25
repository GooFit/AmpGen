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

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace AmpGen;
CoherentSum::CoherentSum() = default; 

CoherentSum::CoherentSum( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix )
  : m_protoAmplitudes( mps )
  , m_evtType( type )
  , m_printFreq(NamedParameter<size_t>(     "CoherentSum::PrintFrequency", 100)  )
  , m_dbThis   (NamedParameter<bool>(       "CoherentSum::Debug"         , false))
  , m_verbosity(NamedParameter<bool>(       "CoherentSum::Verbosity"     , 0)    )
  , m_objCache (NamedParameter<std::string>("CoherentSum::ObjectCache"   ,"")    )
  , m_prefix( prefix )
{
  auto amplitudes      = m_protoAmplitudes.getMatchingRules( m_evtType, prefix);
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
    m_matrixElements[i] = TransitionMatrix<complex_t>( amplitudes[i].first, amplitudes[i].second, mps, this->m_evtType.getEventFormat(), this->m_dbThis);
    CompilerWrapper().compile( m_matrixElements[i].amp, this->m_objCache); } );
  }
  m_isConstant = false;
}

void updateCache(EventList* events, TransitionMatrix<complex_t>& me, const size_t& sizeMax)
{
  if ( me.addressData == 999 )
  { 
    if( events->at(0).cacheSize() <= sizeMax) events->resizeCache(sizeMax);
    me.addressData = events->registerExpression( me.amp );
  }
  events->updateCache(me.amp, me.addressData);
}

void CoherentSum::prepare()
{
  if ( m_weightParam != nullptr ) m_weight = m_weightParam->mean();
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  transferParameters(); 
  std::vector<size_t> changedPdfIndices;
  ProfileClock clockEval; 
  bool print    = false;
  for ( size_t i = 0; i < m_matrixElements.size(); ++i ) {
    m_matrixElements[i].amp.prepare();
    if ( m_prepareCalls != 0 && !m_matrixElements[i].amp.hasExternalsChanged() ) continue;
    if ( m_events != nullptr ) updateCache( m_events, m_matrixElements[i], m_matrixElements.size() ); 
    m_integrator.prepareExpression( m_matrixElements[i].amp );
    changedPdfIndices.push_back(i);
    m_matrixElements[i].amp.resetExternals();
    print = true; 
  }
  clockEval.stop();
  ProfileClock clockIntegral;
  if ( m_integrator.isReady())  updateNorms( changedPdfIndices );
  else if ( m_verbosity ) WARNING( "No simulated sample specified for " << this );
  m_norm = norm();
  if ( m_verbosity && print ) {
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
  //for ( auto& i : changedPdfIndices ) m_integrator.prepareExpression( m_matrixElements[i].amp );
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
          << " A = [ "  << std::real(A)             << " " << std::imag(A) 
          << " ] g = [ "<< std::real(me.coupling()) << " " << std::imag(me.coupling()) << " ]" );
      if( m_dbThis ) me.amp.debug( evt.address() );
      me.coupling.print();
    }
  else
    for ( auto& me : m_matrixElements )
      if ( me.amp.name().find( nameMustContain ) != std::string::npos ) me.amp.debug( evt.address() );
  if( evt.cacheSize() != 0 ) INFO( "Pdf = " << prob_unnormalised( evt ) );
  INFO( "A(x) = " << getVal(evt) );
}

std::vector<FitFraction> CoherentSum::fitFractions(const LinearErrorPropagator& linProp)
{
  prepare();
  bool recomputeIntegrals    = NamedParameter<bool>("CoherentSum::RecomputeIntegrals", false );
  std::vector<FitFraction> outputFractions;
  for(auto& rule : m_protoAmplitudes.rules()) 
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
  auto ffForHead = m_protoAmplitudes.rulesForDecay(m_evtType.mother(), m_prefix);
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
  stream << CompiledExpression< std::complex<double>, const double*, const int&>( amplitude  , "AMP" ) << std::endl; 
  stream << CompiledExpression< double, const double*, const int&>(fcn::norm(amplitude) / normalisation, "FCN" ) << std::endl; 
  if( includePythonBindings ){
    stream << CompiledExpression< unsigned int >( m_matrixElements.size(), "matrix_elements_n" ) << std::endl;
    stream << CompiledExpression< double >      ( normalisation, "normalization") << std::endl;

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
    stream << "extern \"C\" void FCN_mt(double* out, double* events, unsigned int size, int parity, double* amps){\n";
    stream << "  unsigned int n = std::thread::hardware_concurrency();\n";
    stream << "  unsigned int batch_size = size / n;\n";
    stream << "  std::vector<std::thread> threads;\n";
    stream << "  for(size_t i=0; i<n; i++) {\n";
    stream << "    size_t start = batch_size*i;\n";
    stream << "    size_t len = i+1!=n ? batch_size : size-start;\n";
    stream << "    threads.emplace_back(FCN_all, out+start"
      << ", events+start*" 
      << ( *this->m_events )[0].size() 
      << ", len, parity, amps);\n";
    stream << "  }\n";
    stream << "  for(auto &thread : threads)\n";
    stream << "    thread.join();\n";
    stream << "}\n\n";

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
  return std::accumulate( m_matrixElements.begin(), 
      m_matrixElements.end(), 
      complex_t(0,0), 
      [&evt]( auto& a, auto& b ){ return a + b.coefficient * b(evt);} );
}

complex_t CoherentSum::getValNoCache(const Event& evt, const size_t& offset) const
{
  return std::accumulate( m_matrixElements.begin(), 
      m_matrixElements.end(), 
      complex_t(0,0), 
      [&evt,&offset]( auto& a, auto& b ){ return a + b.coefficient * b(evt, offset);} );
}

void CoherentSum::reset( bool resetEvents )
{
  m_prepareCalls                                     = 0;
  m_lastPrint                                        = 0;
  for ( auto& mE : m_matrixElements ) mE.addressData = 999;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integrator = integrator();
  }
}

void CoherentSum::setEvents( EventList& list )
{
  if ( m_verbosity ) INFO( "Setting event list with:" << list.size() << " events for " << this );
  reset();
  m_events = &list;
}

void CoherentSum::setMC( EventList& sim )
{
  if ( m_verbosity ) INFO( "Setting norm. event list with:" << sim.size() << " events for " << this );
  reset();
  m_integrator = integrator(&sim);
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
}

void CoherentSum::printVal(const Event& evt)
{
  for ( auto& mE : m_matrixElements ) {
    unsigned int address = mE.addressData;
    std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << evt.getCache( address )
      << " address = " << address << " " << mE( evt ) << std::endl;
    if( mE.coupling.couplings.size() != 1 ){
      std::cout << "CouplingConstants: " << std::endl;
      mE.coupling.print();
      std::cout << "================================" << std::endl;
    }
  }
}

std::vector<size_t> CoherentSum::cacheAddresses( const EventList& evts ) const
{
  std::vector<size_t> addresses;
  std::transform( m_matrixElements.begin(), m_matrixElements.end(), std::back_inserter(addresses),
      [&evts](auto& it ){ return evts.getCacheIndex( it.amp ) ; } );
  return addresses;
}

complex_t CoherentSum::getVal( const Event& evt ) const
{
  complex_t value( 0., 0. );
  for ( auto& mE : m_matrixElements ) {
    value += mE.coefficient * evt.getCache( mE.addressData );
  }
  return value;
}

complex_t CoherentSum::getVal( const Event& evt, const std::vector<size_t>& cacheAddresses ) const
{
  complex_t value( 0., 0. );
  for ( size_t i = 0; i < m_matrixElements.size(); ++i )
    value += m_matrixElements[i].coefficient * evt.getCache( cacheAddresses[i] );
  return value;
}
