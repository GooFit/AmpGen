#include <Rtypes.h>
#include <TH1.h>
#include <dlfcn.h>
#include <memory>
#include <string>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "TGraph2D.h"



#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#include "AmpGen/DynamicFCN.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/RecursivePhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/QcGenerator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/enum.h"

using namespace AmpGen;

namespace AmpGen { make_enum(generatorType, CoherentSum, PolarisedSum, FixedLib, RGenerator) }

struct FixedLibPDF {
  void* lib = {nullptr};
  AmpGen::DynamicFCN<double( const double*, int )> PDF;

  void prepare(){};
  void setEvents( AmpGen::EventList& evts ){};
  double prob_unnormalised( const AmpGen::Event& evt ) const { return PDF( evt, 1 ); }
  FixedLibPDF( const std::string& lib )
  {
    void* handle = dlopen( lib.c_str(), RTLD_NOW );
    if ( handle == nullptr ) ERROR( dlerror() );
    PDF = AmpGen::DynamicFCN<double( const double*, int )>( handle, "FCN" );
  }
  size_t size() { return 0; }
  void reset( const bool& flag = false ){};
};

template <class PDF_TYPE, class PRIOR_TYPE> 
  void GenerateEvents( EventList& eventsSig
                       , EventList& eventsTag
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& priorSig
                       , PRIOR_TYPE& priorTag
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm )
{
  QcGenerator<PRIOR_TYPE> signalGenerator( priorSig, priorTag );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, eventsSig, eventsTag, nEvents );
}


int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
  auto genType        = NamedParameter<generatorType>( "Type", generatorType::CoherentSum, optionalHelpString("Generator configuration to use:", 
    { {"CoherentSum" , "Full phase-space generator with (pseudo)scalar amplitude"}
    , {"PolarisedSum", "Full phase-space generator with particles carrying spin in the initial/final states"}
    , {"FixedLib"    , "Full phase-space generator with an amplitude from a precompiled library"}
    , {"RGenerator"  , "Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons"} } ) );
  
  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );

  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << nCores  << " / " << concurentThreadsSupported  << " threads" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif
 for( auto& tag : tags ){


    auto tokens       = split(tag, ' ');
    INFO("tag = "<<tokens[0]);
 
  TRandom3 rand;
  rand.SetSeed( seed + 934534 );

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  

  Particle p;
  bool debug=false;
  

//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );






//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );



  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  INFO("Generating time-dependence? " << eventType.isTimeDependent() );
  EventList acceptedSig( sigType );
  EventList acceptedTag( tagType );

  INFO("Generating events with type = " << sigType << " and "<<tagType);


 // if ( genType == generatorType::CoherentSum ) {
    CorrelatedSum cs( sigType, tagType, MPS );
    PhaseSpace phspSig(sigType,&rand);
    PhaseSpace phspTag(tagType,&rand);
    GenerateEvents( acceptedSig, acceptedTag, cs, phspSig, phspTag , nEvents, blockSize, &rand );
    if( acceptedSig.size() == 0 ) return -1;
    TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
    INFO( "Writing output file " );
    acceptedSig.tree(tokens[0])->Write();
    acceptedTag.tree(tokens[0])->Write();


    f->Close();

  //} 
 // else {
 //   FATAL("Did not recognise configuration: " << genType );
 // }
 

 }
  /*
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  acceptedSig.tree( "DalitzEventList" )->Write();
  auto plots = acceptedSig.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(nBins);
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
        acceptedSig.makeProjection( Projection2D(proj[i], proj[j]), LineColor(kBlack) )->Write(); 
      }
    }
  } 
  */




}
