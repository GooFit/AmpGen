#include <Rtypes.h>
#include <TH1.h>
#include <dlfcn.h>
#include <memory>
#include <string>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

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
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"

using namespace AmpGen;

struct FixedLibPDF {
  void* lib;
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


template < class PDF_TYPE, class PRIOR_TYPE > 
  void GenerateEvents( EventList& events
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& prior
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm )
{
  Generator<PRIOR_TYPE> signalGenerator( prior );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, events, nEvents );
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" );
  
  std::string gen_type = NamedParameter<std::string>( "Type", "default", 
      "Generator configuration to use: \n" 
      "\033[3m default            \033[0m: Full phase-space generator with scalar amplitude \n"
      "\033[3m PolarisedSum \033[0m: Full phase-space generator with polarised spin 1/2 amplitude\n"
      "\033[3m FixedLib           \033[0m: Full phase-space generator with an amplitude from a precompiled library\n"
      "\033[3m RGenerator         \033[0m: Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons" );
  
  std::string lib = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");

  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif

  TRandom3 rand;
  rand.SetSeed( seed + 934534 );

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  
  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  EventList accepted( eventType );

  INFO("Generating events with type = " << eventType );

  if ( gen_type == "default" ) {
    CoherentSum sig( eventType, MPS, "" );
    PhaseSpace phsp(eventType,&rand);
    GenerateEvents( accepted, sig, phsp , nEvents, blockSize, &rand );
  } 
  else if ( gen_type == "PolarisedSum" ){
    PolarisedSum sig( eventType, MPS );
    //PhaseSpace phsp(eventType,&rand);
    RecursivePhaseSpace phsp( sig.matrixElements()[0].decayTree->quasiStableTree() , eventType, &rand );
    GenerateEvents( accepted, sig, phsp, nEvents, blockSize, &rand );
  }
  else if ( gen_type == "RGenerator" ) {
    CoherentSum sig( eventType, MPS, "" );
    Generator<RecursivePhaseSpace> signalGenerator( sig[0].decayTree->quasiStableTree(), eventType );
    signalGenerator.setRandom( &rand );
    signalGenerator.fillEventList( sig, accepted, nEvents );
  }
  else if ( gen_type == "FixedLib" ) {
    Generator<> signalGenerator( eventType );
    signalGenerator.setRandom( &rand );
    signalGenerator.setBlockSize( blockSize );
    signalGenerator.setNormFlag( false );
    FixedLibPDF pdf( lib  );
    signalGenerator.fillEventList( pdf, accepted, nEvents );
  } 
  else {
    ERROR("Did not recognise configuration: " << gen_type );
  }
  if( accepted.size() == 0 ) return -1;
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  accepted.tree( "DalitzEventList" )->Write();
  auto plots = accepted.makeDefaultPlots( LineColor(kBlack) );
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(100);
    INFO("Making 2D projections...");
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){
      
        accepted.makeProjection( Projection2D(proj[i],proj[j] ) )->Write(); 
      }
    }
  } 
  
  INFO( "Writing output file " );

  f->Close();
}
