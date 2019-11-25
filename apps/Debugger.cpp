#include <TLorentzVector.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Utilities.h"
#include "TRandom3.h"
#include "AmpGen/AddCPConjugate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"


using namespace AmpGen;

void invertParity( Event& event, const size_t& nParticles)
{
  for( size_t i = 0 ; i < nParticles; ++i )
  {
    event[4*i + 0] = -event[4*i+0];
    event[4*i + 1] = -event[4*i+1];
    event[4*i + 2] = -event[4*i+2];
  }
}

template < class FCN > void debug( FCN& sig, EventList& accepted, bool verbose, TRandom3* rndm, MinuitParameterSet& mps ){
  INFO("Debugging: ");
  sig.setEvents( accepted );
  sig.prepare();
  sig.debug( accepted[0] );
  accepted[0].print();
//  if( verbose ) print( accepted[0], sig.matrixElements(), verbose ); 
  invertParity(accepted[0], accepted.eventType().size() );
  accepted[0].print();
  sig.reset();
  sig.prepare();
  sig.debug( accepted[0] );
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  int seed = NamedParameter<int>( "Seed", 156 );
  TRandom3* rndm = new TRandom3( seed );

  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  bool verbose = NamedParameter<bool>("CoherentSum::Debug", 0 ) || 
                 NamedParameter<bool>("PolarisedSum::Debug", 0 );
  INFO("Using verbose mode: " << verbose );
  AmpGen::MinuitParameterSet MPS;
  MPS.loadFromStream();

  if ( NamedParameter<bool>( "conj", false ) == true ) 
  {
    eventType = eventType.conj();
    INFO( eventType );
    AddCPConjugate(MPS);
  }
  INFO( "EventType = " << eventType );
  
  std::string infile = NamedParameter<std::string>("InputFile","");
  EventList accepted = infile == "" ? EventList( eventType ) : EventList( infile, eventType );
  
  std::string input_units = NamedParameter<std::string>("Units","GeV");
  if( input_units == "MeV" && infile != "") accepted.transform([](auto& event){ for( int i = 0;i<16;++i) event[i]/=1000; } );
  if( infile == "" ){
    Event evt = PhaseSpace( eventType, rndm ).makeEvent();
    accepted.push_back(evt);
  }
  accepted[0].print();

  std::string type = NamedParameter<std::string>("Type","CoherentSum");

  if( type == "PolarisedSum")
  {
    PolarisedSum sig( eventType, MPS );  
    sig.setEvents( accepted );
    sig.prepare();
    debug( sig, accepted, verbose, rndm , MPS );
    sig.setMC(accepted);
    INFO("norm = " << sig.norm() );   
  }
  else if( type == "CoherentSum" )
  {
    CoherentSum sig(eventType, MPS);
    debug(sig, accepted, verbose, rndm, MPS);
    // print(accepted[0], sig.matrixElements() , false);
    INFO( "A(x) = " << sig.getValNoCache( accepted[0] ) );
  }
  else {
    ERROR( "Type: " << type << " is not recognised");
  }
}
