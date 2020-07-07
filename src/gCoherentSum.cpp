#include "AmpGen/gCoherentSum.h"

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
#include "AmpGen/CoherentSum.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace AmpGen;
gCoherentSum::gCoherentSum() = default; 

gCoherentSum::gCoherentSum( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix )
  : m_rules(mps)
  , m_evtTypeA        (type)
  , m_evtTypeC        (type.conj(true))
  , m_printFreq      (NamedParameter<size_t>(     "gCoherentSum::PrintFrequency", 100)  )
  , m_dbThis         (NamedParameter<bool>(       "gCoherentSum::Debug"         , false))
  , m_verbosity      (NamedParameter<bool>(       "gCoherentSum::Verbosity"     , 0)    )
  , m_objCache       (NamedParameter<std::string>("gCoherentSum::ObjectCache"   ,"")    )
  , m_prefix         (prefix)
{
  auto amplitudesA      = m_rules.getMatchingRules( m_evtTypeA, prefix);
  if( amplitudesA.size() == 0 ){
    WARNING("The defined amplitudes don't seem to be able to be able to generate eventType: " << type);
  }
  for( auto& amp : amplitudesA ) INFO( amp.first.decayDescriptor() );
  m_matrixElementsA.resize( amplitudesA.size() );
  m_normalisationsA.resize( m_matrixElementsA.size(), m_matrixElementsA.size() ); 
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , std::thread::hardware_concurrency(), "Number of threads to use" );
  ThreadPool tp(nThreads);
  for(size_t i = 0; i < m_matrixElementsA.size(); ++i){
    tp.enqueue( [i,this,&mps,&amplitudesA]{ 
    m_matrixElementsA[i] = TransitionMatrix<complex_t>( amplitudesA[i].first, amplitudesA[i].second, mps, this->m_evtTypeA.getEventFormat(), this->m_dbThis);
    CompilerWrapper().compile( m_matrixElementsA[i].amp, this->m_objCache); } );
  }
  m_isConstant = false;
  m_x = mps["gCoherentSum::x"];
  m_y = mps["gCoherentSum::y"];
}


void updategCache(EventList* events, TransitionMatrix<complex_t>& me, const size_t& sizeMax)
{
  if ( me.addressData == 999 )
  { 
    if( events->at(0).cacheSize() <= sizeMax) events->resizeCache(sizeMax);
    me.addressData = events->registerExpression( me.amp );
  }
  events->updateCache(me.amp, me.addressData);
}

void gCoherentSum::prepare()
{
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  transferParameters(); 
  std::vector<size_t> changedPdfIndices;
  ProfileClock clockEval; 
  bool print    = false;
  for ( size_t i = 0; i < m_matrixElementsA.size(); ++i ) {
    m_matrixElementsA[i].amp.prepare();
    if ( m_prepareCalls != 0 && !m_matrixElementsA[i].amp.hasExternalsChanged() ) continue;
    if ( m_events != nullptr ) updategCache( m_events, m_matrixElementsA[i], m_matrixElementsA.size() ); 
    m_integratorA.prepareExpression( m_matrixElementsA[i].amp );
    changedPdfIndices.push_back(i);
    m_matrixElementsA[i].amp.resetExternals();
    print = true; 
  }
  clockEval.stop();
  ProfileClock clockIntegral;
  if ( m_integratorA.isReady())  updateNorms( changedPdfIndices );
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

void gCoherentSum::updateNorms( const std::vector<size_t>& changedPdfIndices )
{
  //for ( auto& i : changedPdfIndices ) m_integratorA.prepareExpression( m_matrixElementsA[i].amp );
  std::vector<size_t> cacheIndex;
  std::transform( m_matrixElementsA.begin(), m_matrixElementsA.end(), std::back_inserter(cacheIndex), 
    [this](auto& m){ return this->m_integratorA.getCacheIndex( m.amp ) ; } );
  for ( auto& i : changedPdfIndices )
    for ( size_t j = 0; j < size(); ++j )
      m_integratorA.queueIntegral( cacheIndex[i], cacheIndex[j] ,i, j, &m_normalisationsA );
  m_integratorA.flush();
  m_normalisationsA.resetCalculateFlags();
}

void gCoherentSum::debug( const Event& evt, const std::string& nameMustContain )
{
  prepare();
  for ( auto& me : m_matrixElementsA ) me.amp.resetExternals();
  if ( nameMustContain == "" )
    for ( auto& me : m_matrixElementsA ) {
      auto A = me(evt);
      INFO( std::setw(70) << me.decayTree.uniqueString() 
          << " A = [ "  << std::real(A)             << " " << std::imag(A) 
          << " ] g = [ "<< std::real(me.coupling()) << " " << std::imag(me.coupling()) << " ] "
          << me.decayTree.CP() );

      if( m_dbThis ) me.amp.debug( evt.address() );
    }
  else
    for ( auto& me : m_matrixElementsA )
      if ( me.amp.name().find( nameMustContain ) != std::string::npos ) me.amp.debug( evt.address() );
  if( evt.cacheSize() != 0 ) INFO( "Pdf = " << prob_unnormalised( evt ) );
  INFO( "A(x) = " << getVal(evt) );
}

std::vector<FitFraction> gCoherentSum::fitFractions(const LinearErrorPropagator& linProp)
{
  prepare();
  bool recomputeIntegrals    = NamedParameter<bool>("gCoherentSum::RecomputeIntegrals", false );
  std::vector<FitFraction> outputFractions;
  for(auto& rule : m_rules.rules()) 
  {
    FitFractionCalculator<gCoherentSum> pCalc(this, findIndices(m_matrixElementsA, rule.first), recomputeIntegrals);
    for(auto& process : rule.second) 
    {
      if(process.head() == m_evtTypeA.mother() && process.prefix() != m_prefix) continue;
      auto numeratorIndices   = processIndex(m_matrixElementsA,process.name());
      if(numeratorIndices.size() == 0 || numeratorIndices == pCalc.normSet ) continue; 
      pCalc.emplace_back(process.name(), numeratorIndices);
    }
    if( pCalc.calculators.size() == 0 ) continue;  
    auto fractions = pCalc(rule.first, linProp);
    std::transform( fractions.begin(), fractions.end(), std::back_inserter(outputFractions),[](auto& p){ return p;} );
  };
  auto ffForHead = m_rules.rulesForDecay(m_evtTypeA.mother(), m_prefix);
  FitFractionCalculator<gCoherentSum> iCalc(this, findIndices(m_matrixElementsA, m_evtTypeA.mother()), recomputeIntegrals);
  for ( size_t i = 0; i < ffForHead.size(); ++i ) 
  {
    for ( size_t j = i + 1; j < ffForHead.size(); ++j ) 
    {
      iCalc.emplace_back(ffForHead[i].name() + "x" + ffForHead[j].name(), 
            processIndex(m_matrixElementsA, ffForHead[i].name()), 
            processIndex(m_matrixElementsA, ffForHead[j].name()) );
    }
  }
  std::vector<FitFraction> interferenceFractions = iCalc(m_evtTypeA.mother()+"_interference",linProp);
  for(auto& f : interferenceFractions) outputFractions.push_back(f); 
  for(auto& p : outputFractions) INFO(p);
  return outputFractions;
}

void gCoherentSum::generateSourceCode(const std::string& fname, const double& normalisation, bool add_mt)
{
  std::ofstream stream( fname );
  transferParameters();
  stream << std::setprecision( 10 );
  stream << "#include <complex>\n";
  stream << "#include <vector>\n";
  stream << "#include <math.h>\n";
  if ( add_mt ) stream << "#include <thread>\n";
  bool includePythonBindings = NamedParameter<bool>("gCoherentSum::IncludePythonBindings",false);

  if (includePythonBindings) INFO("Including Python bindings");



  for ( auto& p : m_matrixElementsA ){
    stream << p.amp << std::endl;
    p.amp.compileWithParameters( stream );
    if( includePythonBindings ) p.amp.compileDetails( stream );
  }
  Expression event = Parameter("x0",0,true);
  Expression pa    = Parameter("double(x1)",0,true);
  Expression amplitude;
  for( unsigned int i = 0 ; i < size(); ++i ){
    auto& p = m_matrixElementsA[i];
    Expression this_amplitude = p.coupling() * Function( programatic_name( p.amp.name() ) + "_wParams", {event} ); 
    amplitude = amplitude + ( p.decayTree.finalStateParity() == 1 ? 1 : pa ) * this_amplitude; 
  }
  stream << CompiledExpression< std::complex<double>, const double*, const int&>( amplitude  , "AMP" ) << std::endl; 
  stream << CompiledExpression< double, const double*, const int&>(fcn::norm(amplitude) / normalisation, "FCN" ) << std::endl; 
  if( includePythonBindings ){
    stream << CompiledExpression< unsigned int >( m_matrixElementsA.size(), "matrix_elements_n" ) << std::endl;
    stream << CompiledExpression< double >      ( normalisation, "normalization") << std::endl;

    stream << "extern \"C\" const char* matrix_elements(int n) {\n";
    for ( size_t i = 0; i < m_matrixElementsA.size(); i++ ) {
      stream << "  if(n ==" << i << ") return \"" << m_matrixElementsA.at(i).amp.progName() << "\" ;\n";
    }
    stream << "  return 0;\n}\n";
    stream << "extern \"C\" void FCN_all(double* out, double* events, unsigned int size, int parity, double* amps){\n";
    stream << "  double local_amps [] = {\n";
    for ( unsigned int i = 0; i < size(); ++i ) {
      auto& p = m_matrixElementsA[i];
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
      auto& p    = m_matrixElementsA[i];
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
      auto& p    = m_matrixElementsA[i];
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

complex_t gCoherentSum::getValNoCache( const Event& evt ) const
{
  return std::accumulate( m_matrixElementsA.begin(), 
      m_matrixElementsA.end(), 
      complex_t(0,0), 
      [&evt]( auto& a, auto& b ){ return a + b.coefficient * b(evt);} );
}

complex_t gCoherentSum::getValNoCache(const Event& evt, const size_t& offset) const
{
  return std::accumulate( m_matrixElementsA.begin(), 
      m_matrixElementsA.end(), 
      complex_t(0,0), 
      [&evt,&offset]( auto& a, auto& b ){ return a + b.coefficient * b(evt, offset);} );
}

void gCoherentSum::reset( bool resetEvents )
{
  m_prepareCalls                                     = 0;
  m_lastPrint                                        = 0;
  for ( auto& mE : m_matrixElementsA ) mE.addressData = 999;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integratorA = integrator();
  }
}

void gCoherentSum::setEvents( EventList& list )
{
  if ( m_verbosity ) INFO( "Setting event list with:" << list.size() << " events for " << this );
  reset();
  m_events = &list;
}

void gCoherentSum::setMC( EventList& sim )
{
  if ( m_verbosity ) INFO( "Setting norm. event list with:" << sim.size() << " events for " << this );
  reset();
  m_integratorA = integrator(&sim);
}

real_t gCoherentSum::norm() const
{
  return norm(m_normalisationsA);
}

real_t gCoherentSum::norm(const Bilinears& norms) const
{
  complex_t acc(0, 0);
  for ( size_t i = 0; i < size(); ++i ) {
    for ( size_t j = 0; j < size(); ++j ) {
      auto val = norms.get(i, j) 
        * m_matrixElementsA[i].coefficient 
        * std::conj(m_matrixElementsA[j].coefficient);
      acc += val;
    }
  }
  real_t N = acc.real();
  N = N * pow(m_x->mean(), 2) + pow(m_y->mean(), 2);// + pow(N, 0.5) * m_x->mean() * 2;
  return N;
}

complex_t gCoherentSum::norm(const size_t& x, const size_t& y) const
{
  return m_normalisationsA.get(x, y);
}

void gCoherentSum::transferParameters()
{
  for ( auto& mE : m_matrixElementsA ) mE.coefficient = mE.coupling();
  m_weight.update();
}

void gCoherentSum::printVal(const Event& evt)
{
  for ( auto& mE : m_matrixElementsA ) {
    unsigned int address = mE.addressData;
    std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << evt.getCache( address )
      << " address = " << address << " " << mE( evt ) << std::endl;
    if( mE.coupling.size() != 1 ){
      std::cout << "CouplingConstants: " << std::endl;
      mE.coupling.print();
      std::cout << "================================" << std::endl;
    }
  }
}

std::vector<size_t> gCoherentSum::cacheAddresses( const EventList& evts ) const
{
  std::vector<size_t> addresses;
  std::transform( m_matrixElementsA.begin(), m_matrixElementsA.end(), std::back_inserter(addresses),
      [&evts](auto& it ){ return evts.getCacheIndex( it.amp ) ; } );
  return addresses;
}

complex_t gCoherentSum::getVal( const Event& evt ) const
{
  complex_t value( 0., 0. );
  for ( auto& mE : m_matrixElementsA ) {

    value += mE.coefficient * evt.getCache( mE.addressData );
  }
  value = value * complex_t(m_x->mean(), m_y->mean());
  return value;
}

complex_t gCoherentSum::getVal( const Event& evt, const std::vector<size_t>& cacheAddresses ) const
{
  complex_t value( 0., 0. );
  for ( size_t i = 0; i < m_matrixElementsA.size(); ++i )

    value += m_matrixElementsA[i].coefficient * evt.getCache( cacheAddresses[i] );
  return value;
}
