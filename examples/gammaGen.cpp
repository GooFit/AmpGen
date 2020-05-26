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
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/enum.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AddCPConjugate.h"

using namespace AmpGen;

namespace AmpGen { make_enum(generatorType, CoherentSum, PolarisedSum, FixedLib, RGenerator, TreePhaseSpace) }

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

template <class PDF_TYPE, class PRIOR_TYPE> 
  void gammaGenerateEvents( EventList& events
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& prior
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm )
{
  Generator<PRIOR_TYPE> signalGenerator( prior );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  EventList mc( events.eventType() );
  signalGenerator.fillEventListPhaseSpace( mc, blockSize, pdf.size());
  INFO("Have "<<mc.size()<<" integration events");
  pdf.setMC(mc);
  
  pdf.prepare();
  INFO("Prepared pCoherentSum");
  signalGenerator.fillEventList( pdf, events, nEvents , false);
}


int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  bool outputVals     = NamedParameter<bool>        ("outputVals"     , false, "Print Output values" );
 

  std::string SFtype = NamedParameter<std::string> ("SFType", "Psi3770");
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
  std::string ampFile = NamedParameter<std::string>("ampFile"   , "vals.csv" , "Name of printed values" ); 
  auto genType        = NamedParameter<generatorType>( "Type", generatorType::CoherentSum, optionalHelpString("Generator configuration to use:", 
    { {"CoherentSum"     , "Full phase-space generator with (pseudo)scalar amplitude"}
    , {"PolarisedSum"    , "Full phase-space generator with particles carrying spin in the initial/final states"}
    , {"FixedLib"        , "Full phase-space generator with an amplitude from a precompiled library"}
    , {"RGenerator"      , "Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons"}
    , {"TreePhaseSpace"  , "Recursive phase-space generator with generic handling of intermediate states."} } ) );
  
  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );
  auto BTags           = NamedParameter<std::string>("BTagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();

  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << nCores  << " / " << concurentThreadsSupported  << " threads" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif

  TRandom3 rand;
  rand.SetSeed( seed );

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  
  EventType eventType; 
  std::string decay   = NamedParameter<std::string>("Decay","","Single decay written on the command line"); 
  if( decay != "" )
  {
    Particle p(decay);
    eventType = p.eventType();
    MPS.add(p.decayDescriptor()+"_Re", Flag::Fix, 1., 0);
    MPS.add(p.decayDescriptor()+"_Im", Flag::Fix, 0., 0);
  } 
  else eventType = EventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                  NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );
  
  if ( NamedParameter<bool>( "conj", false ) == true ) {
    eventType = eventType.conj();
    INFO( eventType );
    AddCPConjugate(MPS);
  }

  INFO("Generating time-dependence? " << eventType.isTimeDependent() );


  INFO("Generating events with type = " << eventType );




//    sig2.prepare();
//    CoherentSum sig( eventType, MPS );
//    sig.reset();

//    GenerateEvents( accepted, sig, phsp , nEvents, blockSize, &rand );
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  for (auto& BTag : BTags){

INFO("B DecayType = "<<BTag);
  EventList accepted( eventType );
    PhaseSpace phsp(eventType,&rand);
    

    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag,' ')[1];
    int gammaSign = std::stoi(split(BTag,' ')[2]);
    int nEvents_B = nEvents;
    if (split(BTag, ' ').size() > 3){
      nEvents_B = std::stod(split(BTag, ' ')[3]) * nEvents;
    }
    INFO("GammaSign = "<<gammaSign);

    pCoherentSum sig2( eventType, MPS, B_Pref, gammaSign);
    EventList mc( eventType);
      Generator<PhaseSpace> signalGenerator( phsp );
//      CoherentSum sig(eventType, MPS);
//      sig.prepare();

//    signalGenerator.fillEventListPhaseSpace( mc, blockSize, sig2.size());
 //   INFO("Have "<<mc.size()<<" integration events");

//   sig2.setMC(mc);
//   sig2.prepare();
    gammaGenerateEvents(accepted, sig2, phsp, nEvents_B, blockSize, &rand);
    INFO("Finished generating events");
//  if( accepted.size() == 0 ) return -1;

  accepted.tree( B_Name )->Write();
  bool doPlots = NamedParameter<bool>("doPlots", true);
  if(doPlots){
  auto plots = accepted.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(nBins);
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
        std::stringstream ss;
        ss<<B_Name<<"_s"<<i<<j;
        accepted.makeProjection( Projection2D(proj[i], proj[j]), LineColor(kBlack) )->Write(ss.str().c_str()); 
      }
    }
  } 
  }
  
//    GenerateEvents( accepted, sig2, phsp , nEvents, blockSize, &rand );
  if (outputVals){
    std::ofstream out;
    out.open(ampFile.c_str());
    for (size_t i=0; i < accepted.size(); i++){
      auto eventSig = accepted[i]; 
      auto A = sig2.getVal(eventSig);
      INFO("Output = "<<A);
//      out<<ABCD<<"\n";
      auto sig01 = eventSig.s(0,1);
      auto sig02 = eventSig.s(0,2);
      auto sig12 = eventSig.s(1,2);

      out<<sig01<<"\t"<<sig02<<"\t"<<sig12<<"\t"
         <<A.real()<<"\t"<<A.imag()<<"\n";
    }
    out.close();
  }
  }
f->Write();
  f->Close();
  INFO("Wrote events to "<<outfile);

 return 0; 

}
