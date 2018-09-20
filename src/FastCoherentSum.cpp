#include "AmpGen/FastCoherentSum.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
/// STL
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ratio>

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

#ifdef __USE_OPENMP__
#include <omp.h>
#endif

using namespace AmpGen;

FastCoherentSum::FastCoherentSum( const EventType& type, AmpGen::MinuitParameterSet& mps, const std::string& prefix )
  : m_protoAmplitudes( mps )
  , m_events( nullptr )
  , m_sim( nullptr )
  , m_evtType( type )
  , m_weightParam( nullptr )
  , m_prepareCalls( 0 )
  , m_lastPrint( 0 )
  , m_printFreq( 0 )
  , m_weight( 1 )
  , m_prefix( prefix )
    , m_stateIsGood( true )
{
  bool useCartesian = AmpGen::NamedParameter<bool>( "FastCoherentSum::UseCartesian", true );
  m_dbThis          = AmpGen::NamedParameter<bool>( "FastCoherentSum::Debug", false );
  m_printFreq       = AmpGen::NamedParameter<unsigned int>( "FastCoherentSum::PrintFrequency", 100 );
  m_verbosity       = AmpGen::NamedParameter<unsigned int>( "FastCoherentSum::Verbosity", 0 );
  auto amplitudes  = m_protoAmplitudes.getMatchingRules( m_evtType, prefix, useCartesian );
  for( auto& amp : amplitudes ) addMatrixElement( amp, mps );

  for ( auto& p : m_matrixElements ) p.pdf.compile();
  m_isConstant = isFixedPDF(mps);
  m_normalisations.resize( m_matrixElements.size(), m_matrixElements.size() );
}

void FastCoherentSum::addMatrixElement( std::pair<Particle, Coupling>& particleWithCoupling, const MinuitParameterSet& mps )
{

  auto& protoParticle = particleWithCoupling.first;
  auto& coupling      = particleWithCoupling.second;

  if ( !protoParticle.isStateGood() ) {
    ERROR( "Decay tree not configured correctly for " << protoParticle.uniqueString() );
    m_stateIsGood = false;
    return;
  }
  const std::string name = protoParticle.uniqueString();
  DebugSymbols dbExpressions;
  INFO( "Matrix Element = " << name << ", hash = " << FNV1a_hash(name) );
  const Expression expression = protoParticle.getExpression( m_dbThis ? &dbExpressions : nullptr );
  DEBUG( "Got expression for this tree" );
  for ( auto& mE : m_matrixElements ) {
    if ( name == mE.decayTree->uniqueString() ) return;
  }
  m_matrixElements.emplace_back(
      std::make_shared<Particle>( protoParticle ), coupling,
      CompiledExpression<std::complex<double>, const double*, const double*>(
        expression, "p" + std::to_string(FNV1a_hash(name)), m_evtType.getEventFormat(), m_dbThis ? dbExpressions : DebugSymbols() , &mps ) );
}

void FastCoherentSum::prepare()
{
  if ( m_weightParam != nullptr ) m_weight = m_weightParam->mean();
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  preprepare();
  std::vector<unsigned int> changedPdfIndices;
  auto tStartEval = std::chrono::high_resolution_clock::now();
  bool printed    = false;

  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) {
    auto& pdf = m_matrixElements[i].pdf;
    if ( m_prepareCalls != 0 && !pdf.hasExternalsChanged() ) continue;
    auto t_start = std::chrono::high_resolution_clock::now();
    if ( m_events != nullptr ) {
      if ( m_matrixElements[i].addressData == 999 ){
        m_matrixElements[i].addressData = m_events->registerExpression( pdf );
        DEBUG("Registering expression " << i << " = " << m_matrixElements[i].addressData );
      }
      m_events->updateCache( pdf, m_matrixElements[i].addressData );
    } else if ( i == 0 && m_verbosity ) {
      WARNING( "No data events specified for " << this );
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    auto time  = std::chrono::duration<double, std::milli>( t_end - t_start ).count();
    if ( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ) {
      INFO( pdf.name() << " (t = " << time << " ms, nCalls = " << m_prepareCalls << ", events = " << m_events->size() << ")" );
      printed = true;
    }
    changedPdfIndices.push_back( i );
    pdf.resetExternals();
  }
  auto tStartIntegral = std::chrono::high_resolution_clock::now();

  if ( m_sim != nullptr ) {
    updateNorms( changedPdfIndices );
  } else if ( m_verbosity ) {
    WARNING( "No simulated sample specified for " << this );
  }
  m_norm = norm(); 
  if ( m_verbosity && printed ) {
    auto tNow        = std::chrono::high_resolution_clock::now();
    double timeEval  = std::chrono::duration<double, std::milli>( tStartIntegral - tStartEval ).count();
    double timeIntg  = std::chrono::duration<double, std::milli>( tNow - tStartIntegral ).count();
    double timeTotal = std::chrono::duration<double, std::milli>( tNow - tStartEval ).count();
    INFO( "Time Performance: "
        << "Eval = " << timeEval << " ms"
        << ", Integral = " << timeIntg << " ms"
        << ", Total = " << timeTotal << " ms; normalisation = "  << m_norm );
    m_lastPrint = m_prepareCalls;
  }
  m_prepareCalls++;
}

void FastCoherentSum::updateNorms( const std::vector<unsigned int>& changedPdfIndices )
{
  std::vector<bool> integralHasChanged( size() * size(), 0 );
  for ( auto& i : changedPdfIndices ) {
    auto& pdf = m_matrixElements[i].pdf;
    m_integralDispatch.prepareExpression( pdf );
  }

  for ( auto& i : changedPdfIndices ) {
    for ( unsigned int j = 0; j < size(); ++j ) {
      if ( integralHasChanged[i * size() + j] ) continue;
      integralHasChanged[i * size() + j] = true;
      integralHasChanged[j * size() + i] = true;
      m_integralDispatch.addIntegral( m_matrixElements[i].pdf, m_matrixElements[j].pdf,
          [i, j, this]( const complex_t& val ) {
          this->m_normalisations.set( i, j, val );
          this->m_normalisations.set( j, i, std::conj( val ) );
          } );
    }
  }
  m_integralDispatch.flush(); 
  /*
  if( m_prepareCalls == 0 ){
    for( unsigned int i = 0 ; i < size(); ++i ){
      for( unsigned int j = 0 ; j < size(); ++j ){
        INFO( "N["<<i<<", " << j << "] = " << m_normalisations.get(i,j) << " x " << m_matrixElements[i].coefficient << " x " << std::conj( m_matrixElements[j].coefficient ) );
      }
    }
    auto& event = m_integralDispatch.sim->at(0);
    event.print();
    event.printCache();
    INFO( "evt-prop: " << event.weight() << " " << event.genPdf() );
    double N = norm();
    for( unsigned int i = 0 ; i < size(); ++i){
      auto val = m_normalisations.get( i, i ) * m_matrixElements[i].coefficient * std::conj( m_matrixElements[i].coefficient );
      INFO( "Fraction [ " << m_matrixElements[i].decayTree->uniqueString() << " ] = " << val / N );
    }
    double n2 = 0 ;
    for( auto& evt : * m_integralDispatch.sim ){
      n2 += evt.weight() * prob_unnormalised(evt) / evt.genPdf(); 
      if( prob_unnormalised(evt) > 1e6 ){
        WARNING("Event has very high probability: " << prob_unnormalised(evt) );
        evt.print();
      }
    }
    INFO( "Cross-check normalisation: " << N << " -vs- " << n2 / m_integralDispatch.sim->norm() );
  }
  */
}

void FastCoherentSum::debug( const Event& evt, const std::string& nameMustContain )
{
  prepare();
  for ( auto& pdf : m_matrixElements ) pdf.pdf.resetExternals();
  if ( nameMustContain == "" )
    for ( auto& pdf : m_matrixElements ) {
      auto v = pdf(evt);
      INFO( std::setw(90) << pdf.decayTree->uniqueString() << " = " <<  std::real(v) << " " << std::imag(v) << " cached = " << 
          evt.getCache(  pdf.addressData ) 
          );
      if( m_dbThis ) pdf.pdf.debug( evt );
    }
  else
    for ( auto& pdf : m_matrixElements )
      if ( pdf.pdf.name().find( nameMustContain ) != std::string::npos ) pdf.pdf.debug( evt );

  INFO( "Pdf = " << prob_unnormalised( evt ) );
}

std::map<std::string, std::vector<unsigned int>> FastCoherentSum::getGroupedAmplitudes()
{
  auto rules = m_protoAmplitudes.rulesForDecay( m_evtType.mother() );
  std::map<std::string, std::vector<unsigned int>> ruleMapping;
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) {
    std::string parentProcessName = split( m_matrixElements[i].coupling[0].first->name(), '_' )[0];
    ruleMapping[parentProcessName].push_back( i );
  }
  return ruleMapping;
}

std::vector<FitFraction> FastCoherentSum::fitFractions( const LinearErrorPropagator& linProp )
{
  struct processCalculator {
    std::vector<FFCalculator> calculators;
    std::vector<FitFraction> fractions;
    FitFraction sumV;
    std::string name;
    double sum()
    {
      return std::accumulate( calculators.begin(), calculators.end(), 0., []( double a, auto& b ) { return a + b(); } );
    }

    size_t size() const { return calculators.size() + 1; }
  };

  std::vector<FitFraction> outputFractions;
  std::vector<processCalculator> AllCalculators( m_protoAmplitudes.rules().size() );

  size_t counter = 0;
  size_t pos     = 0;
  for ( auto& processes : m_protoAmplitudes.rules() ) {
    auto& pCalc = AllCalculators[counter++];
    pCalc.name  = processes.first;
    for ( auto& process : processes.second ) {
      INFO( pCalc.name << " " << process.name() );
      if ( process.head() == m_evtType.mother() && process.prefix() != m_prefix ) continue;
      std::string parentProcessName = getParentProcess( process.name() );
      if ( parentProcessName == "" ) continue;
      std::string pName, ppName;
      auto pIndex      = processIndex( process.name() );
      auto parentIndex = processIndex( parentProcessName );
      pCalc.calculators.emplace_back( process.name(), this, pIndex, parentIndex );
    }
    pos += pCalc.size();
  }
  for ( auto& calculator : AllCalculators ) {
    INFO( calculator.name << ": " );
    for ( auto& x : calculator.calculators ) {
      INFO( x.name );
    }
  }

  bool hardcore     = NamedParameter<bool>( "Hardcore", false );
  bool interference = NamedParameter<bool>( "Interference", false );
  auto FitFractions = [this, &AllCalculators, &hardcore]() {
    if ( hardcore )
      this->prepare();
    else
      this->transferParameters();
    std::vector<double> rv;
    for ( auto& pCalc : AllCalculators ) {
      for ( auto& calc : pCalc.calculators ) rv.push_back( calc() );
      rv.push_back( pCalc.sum() );
    }
    return rv;
  };

  auto values = FitFractions();
  auto errors = linProp.getVectorError( FitFractions, values.size() );
  counter     = 0;
  for ( auto& pCalc : AllCalculators ) {
    for ( auto& calc : pCalc.calculators ) {
      pCalc.fractions.emplace_back( calc.name, values[counter], errors[counter] );
      counter++;
    }
    std::sort( pCalc.fractions.begin(), pCalc.fractions.end(),
        []( const FitFraction& f1, const FitFraction& f2 ) { return fabs( f1.val() ) > fabs( f2.val() ); } );
    for ( auto& f : pCalc.fractions ) outputFractions.push_back( f );
    pCalc.sumV = FitFraction( "Sum_" + pCalc.name, values[counter], errors[counter] );
    outputFractions.push_back( pCalc.sumV );
    counter++;
  }
  INFO( "Calculating interference fractions" );
  if ( hardcore && interference ) {
    std::vector<FitFraction> interferenceFractions;
    auto ffForHead = m_protoAmplitudes.rulesForDecay( m_evtType.mother() );

    for ( unsigned int i = 0; i < ffForHead.size(); ++i ) {
      auto process_i = ffForHead[i];
      if ( process_i.prefix() != m_prefix ) continue;
      for ( unsigned int j = i + 1; j < ffForHead.size(); ++j ) {
        auto process_j = ffForHead[j];
        if ( process_j.prefix() != m_prefix ) continue;
        std::string parent_process_name = getParentProcess( process_i.name() );
        FFCalculator iCalc( process_i.name() + "x" + process_j.name(), this, processIndex( process_i.name() ),
            processIndex( process_j.name() ), processIndex( parent_process_name ) );
        interferenceFractions.emplace_back( iCalc.name, iCalc(), linProp.getError( [this, &iCalc]() {
              this->prepare();
              return iCalc();
              } ) );
      }
    }
    std::sort( interferenceFractions.begin(), interferenceFractions.end(),
        []( const FitFraction& f1, const FitFraction& f2 ) { return fabs( f1.val() ) > fabs( f2.val() ); } );
    for ( auto& f : interferenceFractions ) outputFractions.push_back( f );
  }
  for ( auto& p : outputFractions ) {
    INFO( std::setw( 100 ) << p.name() << " " << std::setw( 5 ) << round( p.val() * 100, 3 ) << " ± "
        << round( p.err() * 100, 3 ) << " %" );
  }
  return outputFractions;
}

void FastCoherentSum::generateSourceCode( const std::string& fname, const double& normalisation, bool add_mt )
{
  std::ofstream stream( fname );
  transferParameters();
  stream << std::setprecision( 10 );
  stream << "#include <complex>\n";
  stream << "#include <vector>\n";
  if ( add_mt ) stream << "#include <thread>\n";
  bool include_python_bindings = NamedParameter<bool>("IncludePythonBindings",false);
  for ( auto& p : m_matrixElements ){
    INFO( "Streaming: " << p.decayTree->uniqueString() );
    stream << p.pdf << std::endl;
    INFO("Done streaming, incluing parameters...");
    p.pdf.compileWithParameters( stream );
    if( include_python_bindings ) p.pdf.compileDetails( stream );
  }
  Expression event = Parameter("x0",0,true,0);
  Expression pa    = Parameter("double(x1)",0,true,0);
  Expression amplitude;
  for( unsigned int i = 0 ; i < size(); ++i ){
    auto& p = m_matrixElements[i];
    Expression this_amplitude = p.coupling() * Function( p.pdf.name() + "_wParams", {event} ); 
    amplitude = amplitude + ( p.decayTree->finalStateParity() == 1 ? 1 : pa ) * this_amplitude; 
  }
  
  stream << CompiledExpression< double, const double*, int>(fcn::norm(amplitude) / normalisation, "FCN" ) << std::endl; 
  
  if( include_python_bindings ){
    stream << CompiledExpression< unsigned int >( m_matrixElements.size(), "matrix_elements_n" ) << std::endl;
    stream << CompiledExpression< double >      ( normalisation, "normalization") << std::endl;
    
    stream << "extern \"C\" unsigned int matrix_elements(int n) {\n";
    for ( size_t i = 0; i < m_matrixElements.size(); i++ ) {
      stream << "  if(n ==" << i << ") return " << m_matrixElements.at( i ).pdf.hash() << " ;\n";
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
      int parity = p.decayTree->finalStateParity();
      if ( parity == -1 ) stream << "double(parity) * ";
      stream << "std::complex<double>(amps[" << i * 2 << "],amps[" << i * 2 + 1 << "]) * ";
      stream << p.pdf.name()<< "_wParams( E )";
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
    stream << "    threads.emplace_back(FCN_all, out+start*" << ( *this->m_events )[0].size()
      << ", events+start, len, parity, amps);\n";
    stream << "  }\n";
    stream << "  for(auto &thread : threads)\n";
    stream << "    thread.join();\n";
    stream << "}\n\n";

    stream << "extern \"C\" double coefficients( int n, int which, int parity){\n";
    for ( size_t i = 0; i < size(); i++ ) {
      auto& p    = m_matrixElements[i];
      int parity = p.decayTree->finalStateParity();
      stream << "  if(n == " << i << ") return ";
      if ( parity == -1 ) stream << "double(parity) * ";
      stream << "(which==0 ? " << p.coupling().real() << " : " << p.coupling().imag() << ");\n";
    }
    stream << "  return 0;\n}\n";
  }

  INFO("Generating source-code for PDF: " << fname << " include MT symbols? " << add_mt << " normalisation = " << normalisation );
  stream.close();
}

std::vector<unsigned int> FastCoherentSum::processIndex( const std::string& label ) const
{
  std::vector<unsigned int> indices;
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) {
    bool couplingIncludesThis = false;
    for ( auto& reAndIm : m_matrixElements[i].coupling.couplings )
      if ( reAndIm.first->name().find( label ) != std::string::npos ) couplingIncludesThis = true;
    if ( couplingIncludesThis ) indices.push_back( i );
  }
  return indices;
}

std::string FastCoherentSum::getParentProcess( const std::string& label ) const
{
  auto pI = processIndex( label );
  if ( pI.size() == 0 ) return "";
  auto coupling = m_matrixElements[pI[0]].coupling.couplings;
  for ( unsigned int i = 0; i < coupling.size(); ++i ) {
    if ( coupling[i].first->name().find( label ) != std::string::npos ) {
      return i == 0 ? m_evtType.mother() : coupling[i - 1].first->name();
    }
  }
  return "";
}

unsigned int FastCoherentSum::getPdfIndex( const std::string& name ) const
{
  for ( unsigned int i = 0; i < size(); ++i ) {
    if ( m_matrixElements[i].decayTree->uniqueString() == name ) return i;
  }
  ERROR( "Component " << name << " not found" );
  return 999;
}

bool FastCoherentSum::isFixedPDF(const MinuitParameterSet& mps) const
{
  for ( auto& p : m_matrixElements ) {
    for( auto& c : p.coupling.couplings ) 
      if( c.first->iFixInit() == 0 or 
          c.second->iFixInit() == 0 ) return false; 
  }
  return true;
}

void FastCoherentSum::PConjugate()
{
  for ( auto& amp : m_matrixElements ) {
    if ( amp.decayTree->finalStateParity() == -1 ) {
      auto& top_coupling = *amp.coupling.couplings.begin();
      top_coupling.second->setCurrentFitVal( top_coupling.second->mean() + M_PI );
    }
  }
}

std::complex<double> FastCoherentSum::getValNoCache( const Event& evt ) const
{
  std::complex<double> value( 0, 0 );
  for ( auto& mE : m_matrixElements ) {
    value += mE.coefficient * mE( evt );
  }
  return value;
}

void FastCoherentSum::preprepare()
{
  transferParameters();
  bool isReady = true; 
  for ( auto& mE : m_matrixElements ) {
    if( m_prepareCalls == 0 ) isReady &= mE.pdf.isReady();
    mE.pdf.prepare();
  }
}

void FastCoherentSum::reset( bool resetEvents )
{
  m_prepareCalls                                    = 0;
  m_lastPrint                                       = 0;
  for ( auto& mE : m_matrixElements ) mE.addressInt = 999;
  if ( resetEvents ) {
    m_events = nullptr;
    m_sim    = nullptr;
  }
}
void FastCoherentSum::setEvents( EventList& list )
{
  if ( m_verbosity ) INFO( "Setting events to size = " << list.size() << " for " << this );
  reset();
  m_events = &( list );
}
void FastCoherentSum::setMC( EventList& sim )
{
  if ( m_verbosity ) INFO( "Setting MC = " << &sim << " for " << this );
  if ( m_sim == &sim ) return;
  reset();
  m_sim                  = &sim;
  m_integralDispatch.sim = m_sim;
}

double FastCoherentSum::norm() const
{
  std::complex<double> acc( 0, 0 );
  for ( size_t i = 0; i < size(); ++i ) {
    for ( size_t j = 0; j < size(); ++j ) {
      auto val =
        m_normalisations.get( i, j ) * m_matrixElements[i].coefficient * std::conj( m_matrixElements[j].coefficient );
      acc += val;
    }
  }
  return acc.real();
}

double FastCoherentSum::norm( const Bilinears& norms ) const
{
  std::complex<double> acc( 0, 0 );
  for ( size_t i = 0; i < size(); ++i ) {
    for ( size_t j = 0; j < size(); ++j ) {
      auto val = norms.get( i, j ) * m_matrixElements[i].coefficient * std::conj( m_matrixElements[j].coefficient );
      acc += val;
    }
  }
  return acc.real();
}

std::complex<double> FastCoherentSum::norm( const unsigned int& x, const unsigned int& y ) const
{
  return m_normalisations.get( x, y );
}

void FastCoherentSum::transferParameters()
{
  for ( auto& mE : m_matrixElements ) mE.coefficient = mE.coupling();
}

void FastCoherentSum::printVal( const Event& evt, bool isSim )
{
  for ( auto& mE : m_matrixElements ) {
    unsigned int address = mE.addressData;

    std::cout << mE.decayTree->uniqueString() << " = " << mE.coefficient << " x " << evt.getCache( address )
      << " address = " << address << " " << mE( evt ) << std::endl;
    if( mE.coupling.couplings.size() != 1 ){
      std::cout << "Couplings: " << std::endl;
      for ( auto& ci : mE.coupling.couplings ) {
        std::cout << ci.first->name() << " " << ci.first->mean() << " " << ci.second->mean() << std::endl;
      }
      std::cout << "================================" << std::endl;
    }
  }
}
