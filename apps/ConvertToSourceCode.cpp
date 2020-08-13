#include <complex>
#include <fstream>
#include <string>
#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ThreadPool.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/simd/utils.h"
#include "TRandom3.h"

using namespace AmpGen;

/*
template <class T>
void create_integration_tests(T& pdf, 
    const EventType& type,
    const MinuitParameterSet& mps, 
    const std::vector<Event>& testEvents, 
    const std::string& sourceFile)
{
  auto stringify = [](const std::string arg ){ return "\"" + arg + "\"" ; };  
  std::ofstream unit_tests; 
  unit_tests.open(split( sourceFile, '.')[0] + "_test.cpp");
  unit_tests << "#define BOOST_TEST_DYN_LINK" << std::endl; 
  unit_tests << "#define BOOST_TEST_MODULE amp" << std::endl; 
  unit_tests << "#include <boost/test/included/unit_test.hpp>" << std::endl; 
  unit_tests << "#include \"AmpGen/Particle.h\"" << std::endl; 
  unit_tests << "#include \"AmpGen/CompiledExpression.h\"" << std::endl; 
  unit_tests << "#include \"AmpGen/EventType.h\"" << std::endl; 
  unit_tests << "#include \"AmpGen/MinuitParameterSet.h\"" << std::endl; 
  unit_tests << "#include \"AmpGen/OptionsParser.h\"" << std::endl; 
  unit_tests << "using namespace AmpGen;" << std::endl; 
  
  unit_tests << "void setupOptions(){" << std::endl; 
  for( auto& p : *OptionsParser::getMe() )
  {
    unit_tests << " OptionsParser::setArg( \"" << vectorToString(p.second," ") <<"\");"<< std::endl;    
  }
  unit_tests << "\n}\n" << std::endl; 
 
  for( auto& mE : pdf.matrixElements() ){
    auto value = mE.amp(testEvents[0].address()); 
    unit_tests << "BOOST_AUTO_TEST_CASE( " << mE.amp.progName() + "_test){" << std::endl;
    unit_tests << "  EventType type({" << stringify(type.mother()) << ", " << vectorToString( type.finalStates(), ", ", stringify )  << "});" << std::endl; 
    unit_tests << "  Particle p("<<stringify(mE.decayDescriptor()) << ", type.finalStates());" << std::endl; 
    unit_tests << "  setupOptions();" << std::endl; 
    unit_tests << "  MinuitParameterSet mps; mps.loadFromStream();" << std::endl; 

    unit_tests << "  double event[] = {"<< std::setprecision(15) ;
    for( size_t i = 0 ; i < testEvents[0].size() -1; ++i) unit_tests << testEvents[0][i] << ", ";
    unit_tests << testEvents[0][testEvents[0].size()-1];
    unit_tests << "};" << std::endl; 
    unit_tests << "  auto expr = make_expression<complex_t>(p.getExpression(), p.decayDescriptor(), type.getEventFormat(), mps);" << std::endl; 
    unit_tests << "  auto eval = expr(event);" << std::endl;
    unit_tests << "  BOOST_TEST( std::real(eval) == " << std::real(utils::get<0>(value))<< ", boost::test_tools::tolerance(1e-6)) ;" << std::endl;
    unit_tests << "  BOOST_TEST( std::imag(eval) == " << std::imag(utils::get<0>(value))<< ", boost::test_tools::tolerance(1e-6)) ;" << std::endl;
    unit_tests << "}\n\n";
  }
  unit_tests.close();
}
*/

template <class T> void generate_source(T& pdf, const std::string& sourceFile, MinuitParameterSet& mps, const double& sf)
{
  bool normalise    = NamedParameter<bool>("Normalise",true);
  double safetyFactor = NamedParameter<double>( "SafetyFactor", 3 );
  int seed            = NamedParameter<int>("Seed", 1);
  size_t nEvents      = NamedParameter<size_t>( "NormEvents", 1000000 );
  auto oEventType     = NamedParameter<std::string>("EventType").getVector();
  
  TRandom3 rnd(seed);

  EventType eventType( oEventType ); 
  Generator<PhaseSpace> phsp(eventType);
  phsp.setRandom(&rnd);
  EventList normEvents = phsp.generate(nEvents);
  if constexpr( std::is_same<T, CoherentSum>::value ) pdf.prepare();

  double norm = 1; 
  if( normalise ){
    double pMax = 0;
    for ( auto& evt : normEvents ) 
    {
      if constexpr ( std::is_same<T, PolarisedSum>::value )
      {
        double px, py, pz; 
        rnd.Sphere(px,py,pz, rnd.Uniform(0,1));
        mps["Px"]->setCurrentFitVal(px);
        mps["Py"]->setCurrentFitVal(py);
        mps["Pz"]->setCurrentFitVal(pz);
        pdf.transferParameters();
      }
      double n = 0;
      if constexpr ( std::is_same<T, CoherentSum>::value ) n = std::norm( pdf.getValNoCache(evt) );
      if constexpr ( std::is_same<T, PolarisedSum>::value ) n = pdf.getValNoCache(evt);
      if ( n > pMax ) pMax = n;
    }
    norm = pMax * sf ; 
    INFO( "Making binary with " << pMax << " x safety factor = " << sf );
  }
  mps.resetToInit(); 
  pdf.generateSourceCode( sourceFile, norm, true );
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );
  std::vector<std::string> oEventType = NamedParameter<std::string>( "EventType" ).getVector();
  std::string sourceFile              = NamedParameter<std::string>( "Output" , "output.cpp" );
  std::string type                    = NamedParameter<std::string>( "Type", "CoherentSum" );
  std::string outputPS                = NamedParameter<std::string>( "OutputEvents", "" );
  unsigned int NormEvents             = NamedParameter<unsigned int>( "NormEvents", 1000000 );
  double safetyFactor                 = NamedParameter<double>( "SafefyFactor", 3 );
  unsigned seed                       = NamedParameter<unsigned>("Seed", 0);
  EventType eventType( oEventType );

  AmpGen::MinuitParameterSet MPS; //
  MPS.loadFromStream();
  
  if ( NamedParameter<bool>( "conj", false ) == true ) {
    eventType = eventType.conj(true);
    AddCPConjugate(MPS);
  }
  Generator<PhaseSpace> phsp( eventType );
  TRandom3 rnd;
  rnd.SetSeed( seed  );
  gRandom = &rnd;
  phsp.setRandom( &rnd );

  if( type == "CoherentSum" ){
    CoherentSum sig( eventType, MPS, "" );
    generate_source( sig, sourceFile, MPS, safetyFactor );
  }
  if( type == "PolarisedSum" ){
    PolarisedSum sig( eventType, MPS );
    generate_source( sig, sourceFile, MPS, safetyFactor );
  }
}
