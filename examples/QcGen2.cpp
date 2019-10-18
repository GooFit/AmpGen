#include "AmpGen/CoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/MinuitParameterSet.h"

#include <TFile.h>
#include <TRandom3.h>
#include <TRandom.h>
#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

using namespace AmpGen;
using namespace std::complex_literals;
int main( int argc, char* argv[]){

 OptionsParser::setArgs( argc, argv );
  //OptionsParser::setArgs( argc, argv, "Toy simulation for Quantum Correlated Ψ(3770) decays");
  /* */
  //auto time_wall = std::chrono::high_resolution_clock::now();
  //auto time      = std::clock();
  size_t hwt = std::thread::hardware_concurrency();
  size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
  //double luminosity   = NamedParameter<double>("Luminosity"  , 818.3       , "Luminosity to generate. Defaults to CLEO-c integrated luminosity.");
  size_t nEvents      = NamedParameter<size_t>("nEvents"     , 0           , "Can also generate a fixed number of events per tag, if unspecified use the CLEO-c integrated luminosity.");
  size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
  //bool   poissonYield = NamedParameter<bool  >("PoissonYield", true        , "Flag to include Poisson fluctuations in expected yields (only if nEvents is not specified)");
  //bool   noQCFit      = NamedParameter<bool  >("noQCFit"     , false       , "Treat Signal and Tag as uncorrelated and fit the data individually");
  double crossSection = NamedParameter<double>("CrossSection", 3.26 * 1000 , "Cross section for e⁺e⁻ → Ψ(3770) → DD'");
  std::string output  = NamedParameter<std::string>("Output" , "ToyMC.root", "File containing output events"); 
  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
  bool m_debug        = NamedParameter<bool>("Debug", false, "Debug QcFitter output");
    bool doDebugNorm  = NamedParameter<bool>("doDebugNorm", false, "Debug the normalisation of the pdf");
    int nBins = NamedParameter<int>("nBins", 100, "number of bins for projection");
    int nFits = NamedParameter<int>("nFits", 1, "number of repeats of mini.doFits() for debug purposes!");


  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");

for( auto& tag : tags )
 {

       MinuitParameterSet MPS;
  MPS.loadFromStream();
    EventType signalType( pNames );
    auto tokens       = split(tag, ' ');
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType(); 
  
    EventList eventsSig = Generator<>(signalType, &rndm).generate(2e6); 
    EventList eventsTag = Generator<>(tagType, &rndm).generate(2e6); 


    
    CorrelatedSum cs(signalType, tagType, MPS);
    cs.SetEvents(eventsSig, eventsTag);
    PhaseSpace sigphsp = PhaseSpace(signalType);
    PhaseSpace tagphsp = PhaseSpace(tagType);

    for (int i=0; i < nEvents; i++){
        auto sigEvent = sigphsp.makeEvent();
        auto tagEvent = tagphsp.makeEvent();
        //INFO("cs = "<<cs(sigEvent, tagEvent));
    }
    //auto events = generatePHSP(10);
 }

    return 0;
}