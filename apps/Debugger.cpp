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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AmpGen/EventType.h"
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedAmplitude.h"

using namespace AmpGen;

int invert_parameter( AmpGen::MinuitParameter* param, MinuitParameterSet& mps )
{
  const std::string name = param->name();
  size_t pos             = 0;
  std::string prefix     = "";
  std::string new_name   = name;
  int sgn                = 1;
  if ( name.find( "::" ) != std::string::npos ) {
    pos                              = name.find( "::" );
    auto props                       = AmpGen::ParticlePropertiesList::get( name.substr( 0, pos ), true );
    if ( props != nullptr ) new_name = props->anti().name() + name.substr( pos );
  } else {
    auto tokens        = split( name, '_' );
    std::string reOrIm = *tokens.rbegin();
    std::string prefix = "";
    std::string name   = tokens[0];
    if ( tokens.size() == 3 ) {
      prefix = tokens[0];
      name   = tokens[1];
    }
    if ( reOrIm == "Re" || reOrIm == "Im" ) {
      std::vector<std::string> final_state_ordering;
      AmpGen::Particle test( name, final_state_ordering );
      sgn                       = test.conjugate( true );
      if ( reOrIm == "Re" ) sgn = 1;
      new_name                  = ( prefix != "" ? prefix + "_" : "" ) + test.uniqueString() + "_" + reOrIm;
    } else if ( tokens.size() == 2 ) {
      auto props                       = AmpGen::ParticlePropertiesList::get( name );
      if ( props != nullptr ) new_name = props->anti().name() + "_" + tokens[1];
    }
  }
  DEBUG( param->name() << " →  " << new_name << " sgn = " << sgn );
  mps.map().erase( param->name() );
  param->setName( new_name );
  mps.map().emplace( param->name(), param );
  if ( sgn == -1 ) param->setCurrentFitVal( param->mean() + M_PI );
  return sgn;
}

template <class MatrixElements> 
void print( const Event& event, const MatrixElements& matrixElements, bool verbose )
{
  for ( auto& mE : matrixElements ) {
    INFO( mE.decayTree->uniqueString() << " " << mE.coupling() );
    auto terms = mE.coupling.couplings;
    if ( verbose ) {
      for ( auto& term : terms ) {
        INFO( "--> " << term.first->name() << " = (" << term.first->mean() * cos( term.second->mean() ) << " + i " << term.first->mean() * sin( term.second->mean() ) << ")" );
      }
      mE.pdf.debug( event );
    }
  }
}

template < class FCN > void debug( FCN& sig, EventList& accepted, bool verbose, TRandom3* rndm, MinuitParameterSet& mps ){
  INFO("Debugging: ");
  sig.setEvents( accepted );
  sig.prepare();
  sig.debug( accepted[0] );
  Event accepted_boosted = accepted[0];
  boost( accepted_boosted,  { rndm->Uniform(-1,1), rndm->Uniform(-1,1), rndm->Uniform(-1,1)}, 0.99 );
  
  print( accepted[0], sig.matrixElements(), verbose );
  accepted[0].print();
  INFO( "A(x)  = " << sig.prob_unnormalised( accepted[0] ) );
  
  auto unboosted_value = sig.getValNoCache( accepted[0] );
  INFO("Boosted: ");
  accepted_boosted.print();
  auto boosted_value = sig.getValNoCache( accepted_boosted );
  INFO( "A(x)  = " << unboosted_value ) ;
  INFO( "A(Λx) = " << boosted_value );   
  
  accepted[0].invertParity(4);
  accepted[0].print();
  INFO( "A(Px) = " << sig.getValNoCache( accepted[0] ) );
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  int seed = NamedParameter<int>( "Seed", 156 );

  EventType eventType( NamedParameter<std::string>( "EventType" ).getVector() );

  bool verbose = NamedParameter<bool>( "FastCoherentSum::Debug", 0 );
  INFO("Using verbose mode: " << verbose );
  AmpGen::MinuitParameterSet MPS;
  MPS.loadFromStream();
  if ( NamedParameter<bool>( "conj", false ) == true ) {
    eventType = eventType.conj( false );
    for ( auto& param : MPS ) invert_parameter( param, MPS );
  }
  INFO( "EventType = " << eventType );
  
  std::string infile = NamedParameter<std::string>("InputFile","");
  EventList accepted = infile == "" ? EventList( eventType ) : EventList( infile, eventType );

  TRandom3* rndm = new TRandom3( seed );
  if( infile == "" ){
    Event evt = PhaseSpace( eventType, rndm ).makeEvent();
    accepted.push_back(evt);
  }
  accepted[0].print();

  std::string type = NamedParameter<std::string>("Type","FastCoherentSum");

  if( type == "PolarisedAmplitude"){
    PolarisedAmplitude sig( eventType, MPS );  
    sig.setEvents( accepted );
    sig.prepare();
    debug( sig, accepted, verbose, rndm , MPS );
    sig.setMC(accepted);
    INFO("norm = " << sig.norm() );   
  }
  if( type == "FastCoherentSum" ){
    FastCoherentSum sig( eventType, MPS );  
    debug( sig, accepted, verbose, rndm , MPS );
    print( accepted[0], sig.matrixElements() , false ); 
  }
}
