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

#include "TRandom3.h"

using namespace AmpGen;

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
    unit_tests << "  BOOST_TEST( std::real(eval) == " << std::real(value)<< ", boost::test_tools::tolerance(1e-6)) ;" << std::endl;
    unit_tests << "  BOOST_TEST( std::imag(eval) == " << std::imag(value)<< ", boost::test_tools::tolerance(1e-6)) ;" << std::endl;
    unit_tests << "}\n\n";
  }
  unit_tests.close();
}

void invert( MinuitParameter* param, MinuitParameterSet& mps )
{
  const std::string name = param->name();
  size_t pos             = 0;
  std::string new_name   = name; 
  int         sgn        = 1;
  std::string cartOrPolar = NamedParameter<std::string>("CouplingConstant::Coordinates" ,"cartesian");

  if( name.find("::") != std::string::npos ){
    pos = name.find("::");
    auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
    if( props != 0 ) new_name = props->anti().name() + name.substr(pos); 
  }
  else { 
    auto tokens=split(name,'_');
    std::string reOrIm = *tokens.rbegin();
    std::string name   = tokens[0];
    if ( reOrIm == "Re" || reOrIm == "Im" ){
      Particle test = Particle(name).conj();
      if( cartOrPolar == "polar" )     sgn = reOrIm == "Re" ? test.quasiCP() : 1; 
      if( cartOrPolar == "cartesian" ) sgn = test.quasiCP();
      new_name = test.uniqueString() +"_"+reOrIm;
    }
    else if( tokens.size() == 2 ) {
      auto props = AmpGen::ParticlePropertiesList::get( name );
      if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
    }
  }
  mps.rename( param->name(), new_name );
  if( sgn == -1 ){ param->setInit( -1*param->mean() ) ; param->setCurrentFitVal( -1 * param->mean() );}
}

template <class T> void generate_source(T& pdf, EventList& normEvents, const std::string& sourceFile, MinuitParameterSet& mps, const double& sf)
{
  bool normalise                      = NamedParameter<bool>("Normalise",true);
  std::string type                    = NamedParameter<std::string>( "Type", "CoherentSum" );

  double norm = 1; 
  if( normalise ){
    double pMax = 0 ;
    pdf.setEvents( normEvents );
    pdf.prepare();
    pdf.debug( normEvents[0] );
    for ( auto& evt : normEvents ) {
      if( type == "PolarisedSum" ){ 
        double px, py, pz; 
        gRandom->Sphere(px,py,pz, gRandom->Uniform(0,1));
        mps["Px"]->setCurrentFitVal(px);
        mps["Py"]->setCurrentFitVal(py);
        mps["Pz"]->setCurrentFitVal(pz);
        pdf.transferParameters();
      }
      double n = pdf.prob_unnormalised( evt );
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

  EventType eventType( oEventType );

  AmpGen::MinuitParameterSet MPS; //
  MPS.loadFromStream();
  
  if ( NamedParameter<bool>( "conj", false ) == true ) {
    eventType = eventType.conj(true);
    for ( auto& param : MPS ) invert( param, MPS );
  }
  Generator<PhaseSpace> phsp( eventType );
  TRandom3 rnd;

  gRandom = &rnd;
  phsp.setRandom( &rnd );

  EventList phspEvents( oEventType );
  phsp.fillEventListPhaseSpace( phspEvents, NormEvents );

  if( type == "CoherentSum" ){
    CoherentSum sig( eventType, MPS, "" );
    generate_source( sig, phspEvents, sourceFile, MPS, safetyFactor );
    create_integration_tests(sig, eventType, MPS, {phspEvents[15]}, sourceFile );
  }
  if( type == "PolarisedSum" ){
    PolarisedSum sig( eventType, MPS );
    generate_source( sig, phspEvents, sourceFile, MPS, safetyFactor );
  }
  if ( outputPS != "" ) {
    std::ofstream ofile( outputPS );
    ofile << "0x,0y,0z,0t,1x,1y,1z,1t,2x,2y,2z,2t,3x,3y,3z,3t\n";
    for ( auto& event : phspEvents ) {
      for ( size_t i = 0; i < event.size(); i++ ) ofile << ( i == 0 ? "" : "," ) << event[i];
      ofile << "\n";
    }
  }
}
