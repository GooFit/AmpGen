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
#include <nlohmann/json.hpp>

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
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"


using namespace AmpGen;

void invertParity( Event& event, const size_t& nParticles) {
  for( size_t i = 0 ; i < nParticles; ++i )
  {
    event[4*i + 0] = -event[4*i+0];
    event[4*i + 1] = -event[4*i+1];
    event[4*i + 2] = -event[4*i+2];
  }
}

void randomBoost( Event& event, TRandom3* rndm ) {
  auto v = std::make_tuple( rndm->Uniform(), rndm->Uniform(), rndm->Uniform() );
  auto beta = rndm->Uniform();
  Event pol(4); 
  pol[0] =  -0.104734; 
  pol[1] =   -0.328773; 
  pol[2] = 0.617203; 
  pol[3] = 0.0; 
  boost( pol, v, beta ); 
  pol.print();
  boost( event, v, beta); 
}

void randomRotation( Event& event, TRandom3* rndm ) {
  rotate( event, std::make_tuple( rndm->Uniform(), rndm->Uniform(), rndm->Uniform() ), rndm->Uniform() ); 
}

template <typename fcn_type> void writeRefFile ( const std::string& filename, fcn_type& fcn, const EventList& events ){
  using json = nlohmann::json;
  json output;   
  for( unsigned i = 0 ; i != events.size(); ++i ){
    output["event_"+std::to_string(i)] = events[i].data(); 
    output["pdf_"+std::to_string(i)]   = fcn( events[i] );
    for( const auto& elem : fcn.matrixElements() ){
      auto indices = fcn.cache().find( elem.name() );
      std::vector<double> this_cache; 
      for( const auto& index : indices ){
        auto v = utils::at( fcn.cache()( i / utils::size<real_v>::value, index ) , i % utils::size<real_v>::value );
        this_cache.emplace_back( std::real(v) );
        this_cache.emplace_back( std::imag(v) );
      }
      output[elem.name() + "_"+std::to_string(i)] = this_cache; 
    }
  }
  std::ofstream os(filename); 
  os << output.dump(4) << std::endl; 
  os.close();  
}

template <typename FCN> void debug( FCN& sig, EventList& accepted){
  INFO("Debugging: ");
  unsigned eventToDebug = 2;
  sig.setEvents( accepted );
  sig.prepare();
  accepted[eventToDebug].print();
  sig.debug( accepted[eventToDebug] );
  
  /*
  INFO("Parity: " ); 

  for( unsigned int i = 0 ; i != accepted.size(); ++i ) 
    invertParity(accepted[i], accepted.eventType().size() );
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  */
  
  INFO("Random boost:" );

    randomBoost( accepted[eventToDebug],  new TRandom3(eventToDebug) ); 
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  
  INFO("Random rotation:" ); 

  randomRotation( accepted[eventToDebug],  new TRandom3(eventToDebug) ); 
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  
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
  if( NamedParameter<bool>( "AddConj", false) == true && NamedParameter<bool>( "conj", false ) == false )
  {
    AddCPConjugate(MPS); 
  }
  INFO( "EventType = " << eventType );
  
  std::string infile = NamedParameter<std::string>("InputFile","");
  EventList accepted = infile == "" ? EventList( eventType ) : EventList( infile, eventType );
  std::string refFileOutput = NamedParameter<std::string>("RefFileOutput",""); 

  std::string input_units = NamedParameter<std::string>("Units","GeV");
  if( input_units == "MeV" && infile != "") accepted.transform([](auto& event){ for( unsigned i = 0;i< event.size();++i) event[i]/=1000; } );
  if( infile == "" ){
    accepted = Generator<PhaseSpace>(eventType, rndm).generate(16);
    for( unsigned i = 0 ; i != 16; ++i ) accepted[i].setIndex(i);
  }
  std::vector<double> event = NamedParameter<double>("Event",0).getVector();
  if( event.size() != 1 ) accepted[0].set( event.data() );
  
  std::string type = NamedParameter<std::string>("Type","CoherentSum");

  if( type == "PolarisedSum")
  {
    PolarisedSum sig( eventType, MPS );  
    sig.setEvents( accepted );
    sig.prepare();
    if( refFileOutput != "" ) writeRefFile( refFileOutput, sig, accepted ); 
    debug( sig, accepted);
    sig.setMC(accepted);
    INFO("norm = " << sig.norm() );   
  }
  else if( type == "CoherentSum" )
  {
    CoherentSum sig(eventType, MPS);
    sig.setEvents(accepted); 
    sig.prepare(); 
    if( refFileOutput != "" ) writeRefFile( refFileOutput, sig, accepted ); 
    debug(sig, accepted);
    INFO( "A(x) = " << sig.getValNoCache( accepted[0] ) );
  }
  else if( type == "IncoherentSum" )
  {
    IncoherentSum sig(eventType, MPS,"Inco"); 
    sig.setMC(accepted);
    sig.prepare();
    debug(sig, accepted);  
    INFO("norm = " << sig.norm() );   
  }
  else {
    ERROR( "Type: " << type << " is not recognised");
  }
}
