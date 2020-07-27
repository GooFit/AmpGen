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
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/enum.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AddCPConjugate.h"

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_t = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_t = AmpGen::EventList; 
#endif

using namespace AmpGen;

namespace AmpGen { 
  make_enum(pdfTypes, CoherentSum, PolarisedSum, FixedLib)
  make_enum(phspTypes, PhaseSpace, RecursivePhaseSpace, TreePhaseSpace)
}  

struct FixedLibPDF 
{
  void* lib = {nullptr};
  DynamicFCN<double( const double*, int )> PDF;
  void debug( const Event& event) {};
  void prepare(){};
  void setEvents( AmpGen::EventList& evts ){};
  void setEvents( AmpGen::EventListSIMD& evts ){};
  double operator()( const AmpGen::Event& evt ) const { return PDF( evt, 1 ); }
  double operator()( const double* evt, const unsigned& index )
  {
    return PDF(evt, 1 ); 
  } 
  FixedLibPDF( const std::string& lib )
  {
    void* handle = dlopen( lib.c_str(), RTLD_NOW );
    if ( handle == nullptr ) ERROR( dlerror() );
    PDF = DynamicFCN<double( const double*, int )>( handle, "FCN" );
  }
  size_t size() { return 0; }
  void reset( const bool& flag = false ){};
};


template <typename pdf_t> Particle getTopology(const pdf_t& pdf)
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    FATAL("Cannot deduce decay topology from a compiled library, check generator options");
  }
  else return pdf.matrixElements()[0].decayTree.quasiStableTree();
}

template <typename pdf_t> std::vector<Particle> getDecayChains( const pdf_t& pdf )
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    FATAL("Cannot deduce decay topology from a compiled library, check generator options");
  }
  else {
    std::vector<Particle> channels; 
    for( auto& chain : pdf.matrixElements() ) channels.push_back( chain.decayTree );
    return channels;
  }
}

template <typename pdf_t> void generateEvents( EventList& events
                       , pdf_t& pdf 
                       , const phspTypes& phsp_type
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm
                       , const bool& normalise = true )
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    Generator<PhaseSpace, EventList> signalGenerator(events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else { 
  if( phsp_type == phspTypes::PhaseSpace )
  {
    Generator<PhaseSpace, EventList_t> signalGenerator(events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else if( phsp_type == phspTypes::RecursivePhaseSpace )
  {
    Generator<RecursivePhaseSpace, EventList_t> signalGenerator( getTopology(pdf), events.eventType(), rndm );
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents);
  }
  else if( phsp_type == phspTypes::TreePhaseSpace )
  {
    Generator<TreePhaseSpace, EventList_t> signalGenerator(getDecayChains(pdf), events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else {
    FATAL("Phase space configuration: " << phsp_type << " is not supported");
  }
  }
}



int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 5000000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
  auto pdfType        = NamedParameter<pdfTypes>( "Type", pdfTypes::CoherentSum, optionalHelpString("Type of PDF to use:", 
    { {"CoherentSum"     , "Describes decays of a (pseudo)scalar particle to N pseudoscalars"}
    , {"PolarisedSum"    , "Describes the decay of a particle with spin to N particles carrying spin."}
    , {"FixedLib"        , "PDF to describe a decay from a precompiled library, such as those provided to GAUSS."}} ) ); 
  auto phspType        = NamedParameter<phspTypes>( "PhaseSpace", phspTypes::PhaseSpace, optionalHelpString("Phase-space generator to use:", 
    { {"PhaseSpace"         , "Phase space generation based on Raubold-Lynchi algorithm (recommended)."}
    , {"TreePhaseSpace"     , "Divides the phase-space into a series of quasi two-body phase-spaces for efficiently generating narrow states."}
    , {"RecursivePhaseSpace", "Includes possible quasi-stable particles and the phase spaces of their decay products, such as Λ baryons."}} ) ); 
  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );

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
 
  bool conj = NamedParameter<bool>("Conj",false, "Flag to generate the CP conjugate amplitude under the assumption of CP conservation");
  if ( conj == true ) {
    eventType = eventType.conj();
    INFO( eventType );
    AddCPConjugate(MPS);
  }
  if( OptionsParser::printHelp() ) return 0; 
  
  INFO("Writing output: " << outfile );
  #ifdef _OPENMP
  INFO("Using: " << nCores  << " / " << concurentThreadsSupported  << " threads" );
  #endif

  INFO("Generating time-dependence? " << eventType.isTimeDependent() );
  EventList accepted( eventType );

  INFO("Generating events with type = " << eventType );

  if ( pdfType == pdfTypes::CoherentSum ){
    CoherentSum pdf( eventType, MPS);
    generateEvents(accepted, pdf, phspType , nEvents, blockSize, &rand );
  }
  else if ( pdfType == pdfTypes::PolarisedSum ){
    PolarisedSum pdf(eventType, MPS);
    generateEvents( accepted, pdf, phspType, nEvents, blockSize, &rand );
  }
  else if ( pdfType == pdfTypes::FixedLib ){
    FixedLibPDF pdf(lib);
    generateEvents( accepted, pdf, phspType, nEvents, blockSize, &rand, false );
  }
  else {
    FATAL("Did not recognise configuration: " << pdfType );
  }
  if( accepted.size() == 0 ) return -1;
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  accepted.tree( "DalitzEventList" )->Write();
  auto plots = accepted.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(nBins);
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
        accepted.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack) )->Write(); 
      }
    }
  } 
  INFO( "Writing output file " );

  f->Close();
}
