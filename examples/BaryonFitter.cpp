#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <THStack.h>
#include <TCanvas.h>
#include <TFile.h>

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/ErrorPropagator.h"
#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList;
#endif
#include "AmpGen/EventType.h"
#include "AmpGen/Factory.h"
#include "AmpGen/RecursivePhaseSpace.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/Kinematics.h"

#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

using namespace AmpGen;

void randomizeStartingPoint( MinuitParameterSet& mps, TRandom3& rand)
{
  for (auto& param : mps) {
    if ( ! param->isFree() || param->name() == "Px" || param->name() == "Py" || param->name() == "Pz" ) continue;
    double min = param->minInit();
    double max = param->maxInit();
    double new_value = rand.Uniform(param->mean()-param->stepInit(),param->mean()+param->stepInit());
    if( min != 0 && max != 0 )
      new_value = rand.Uniform(min,max);
    param->setInit( new_value );
    param->setCurrentFitVal( new_value );
    INFO( param->name() << "  = " << param->mean() << " " << param->stepInit() );
  }
}

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList_type& data, EventList_type& mc, MinuitParameterSet& MPS )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();

  pdf.setEvents( data );

  /* Minimiser is a general interface to Minuit1/Minuit2,
     that is constructed from an object that defines an operator() that returns a double
     (i.e. the likielihood, and a set of MinuitParameters. */
  Minimiser mini( pdf, &MPS );
  mini.doFit();
  FitResult* fr = new FitResult(mini);

  /* Estimate the chi2 using an adaptive / decision tree based binning,
     down to a minimum bin population of 15, and add it to the output.*/
  //if(data.eventType().size() < 5){
  //  Chi2Estimator chi2( data, mc, pdf, 15 );
  //  //chi2.writeBinningToFile("chi2_binning.txt");
  //  fr->addChi2( chi2.chi2(), chi2.nBins() );
  //}

  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );
  fr->print();

  /* Save weighted data and norm MC for the different components in the PDF, i.e. the signal and backgrounds.
     The structure assumed the PDF is some SumPDF<T1,T2,...>. */
  unsigned int counter = 1;
  for_each(pdf.pdfs(), [&]( auto& f ){
    mc.transform([&f](auto& mcevt){mcevt.setWeight(f.getValNoCache(mcevt)*mcevt.weight()/mcevt.genPdf());}).tree(counter>1?"MCt"+std::to_string(counter):"MCt")->Write();
    data.tree(counter>1?"t"+std::to_string(counter):"t")->Write();
    counter++;
  } );

  return fr;
}

void invertParity( Event& event, const size_t& nParticles=0)
{
  for( size_t i = 0 ; i < nParticles; ++i )
  {
    event[4*i + 0] = -event[4*i+0];
    event[4*i + 1] = -event[4*i+1];
    event[4*i + 2] = -event[4*i+2];
  }
}

int main( int argc, char* argv[] )
{
  gErrorIgnoreLevel = 1001;

  OptionsParser::setArgs( argc, argv );

  const std::vector<std::string> dataFile = NamedParameter<std::string>("DataSample","").getVector();
  const std::string simFile               = NamedParameter<std::string>("SimFile", ""   , "Name of file containing simulated sample for using in MC integration");
  const std::string logFile               = NamedParameter<std::string>("LogFile","Fitter.log");
  const std::string plotFile              = NamedParameter<std::string>("Plots","plots.root");
  const std::string prefix                = NamedParameter<std::string>("PlotPrefix","");
  const std::string idbranch              = NamedParameter<std::string>("IDBranch","");
  const std::string mcidbranch            = NamedParameter<std::string>("MCIDBranch","");
  const std::string weight_branch         = NamedParameter<std::string>("WeightBranch","","Name of branch containing event weights.");
  const std::string mc_weight_branch      = NamedParameter<std::string>("MCWeightBranch","","Name of branch containing event weights.");

  const auto nev_MC = NamedParameter<int>("NEventsMC", 8e6, "Number of MC events for normalization.");
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>(),
                                            "List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto MCbNames = NamedParameter<std::string>("MCBranches", std::vector<std::string>(),
                                              "List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto pNames = NamedParameter<std::string>("EventType" , ""
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector();

#if ENABLE_AVX
  if(!idbranch.empty() || !weight_branch.empty() || !mcidbranch.empty() || !mc_weight_branch.empty()){
    ERROR("Vectorized version currently not supported when adding extra branches");
    return 1;
  }
#endif

#ifdef _OPENMP
  unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
  unsigned int nThreads                  = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported );
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif


  /* A MinuitParameterSet is (unsurprisingly) a set of fit parameters, and can be loaded from
     the parsed options. For historical reasons, this is referred to as loading it from a "Stream" */
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  TRandom3 rndm = TRandom3( NamedParameter<unsigned int>("Seed", 1 ) ) ;
  if( NamedParameter<bool>("RandomizeStartingPoint",false) ) randomizeStartingPoint(MPS,rndm );

  /* An EventType specifies the initial and final state particles as a vector that will be described by the fit.
     It is typically loaded from the interface parameter EventType. */
  EventType evtType(pNames);

  /* A CoherentSum is the typical amplitude to be used, that is some sum over quasi two-body contributions
     weighted by an appropriate complex amplitude. The CoherentSum is generated from the couplings described
     by a set of parameters (in a MinuitParameterSet), and an EventType, which matches these parameters
     to a given final state and a set of data. A common set of rules can be matched to multiple final states,
     i.e. to facilitate the analysis of coupled channels. */
  PolarisedSum sig(evtType, MPS);

  /* Events are read in from ROOT files. If only the filename and the event type are specified,
     the file is assumed to be in the specific format that is defined by the event type,
     unless the branches to load are specified in the user options */
  EventList_type events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch), ExtraBranches(std::vector<std::string>{idbranch}) );

  /* Generate events to normalise the PDF with. This can also be loaded from a file,
     which will be the case when efficiency variations are included. */
  EventList_type eventsMC = simFile == ""
   ? EventList_type(Generator<RecursivePhaseSpace, EventList>(sig.matrixElements()[0].decayTree.quasiStableTree(), events.eventType(), &rndm).generate(nev_MC))
   : EventList_type(simFile, evtType, Branches(MCbNames), WeightBranch(mc_weight_branch), ExtraBranches(std::vector<std::string>{mcidbranch}));

  /* Transform data if we have an ID brach. That branch indicates that we operate on a sample with particles+antiparticles mixed.
     The transformation also includes boosting to the restframe of the head of the decay chain.
     TODO: There might be situations where you want to separate both transformations */
  auto const n_final_state_particles = evtType.size();
  std::vector<unsigned> daughters_as_ints(n_final_state_particles);
  std::iota (daughters_as_ints.begin(), daughters_as_ints.end(), 0u);
  auto frame_transform = [&daughters_as_ints, &n_final_state_particles](auto& event){
    TVector3 pBeam(0,0,1);
    if( event[event.size()-1] < 0 ){
      invertParity( event, n_final_state_particles);
      pBeam = -pBeam;
    }
    TLorentzVector pP = pFromEvent(event,daughters_as_ints);
    //if( pP.P() < 10e-5) return;
    TVector3 pZ = pP.Vect();
    rotateBasis( event, (pBeam.Cross(pZ) ).Cross(pZ), pBeam.Cross(pZ), pZ );
    boost( event, {0, 0, -1}, pP.P()/pP.E() );
  };


  for( auto& event : events )
    if( event[event.size()-1] < 0 ){
      event.print();
      break;
    }
  events.transform( frame_transform );
  for( auto& event : events )
    if( event[event.size()-1] < 0 ){//E if there's no ID branch
      event.print();
      break;
    }
  for( auto& event : eventsMC )
    if( event[event.size()-1] < 0 ){
      event.print();
      break;
    }
  eventsMC.transform( frame_transform );
  for( auto& event : eventsMC )
    if( event[event.size()-1] < 0 ){//E if there's no ID branch
      event.print();
      break;
    }
  sig.setMC(eventsMC);
  sig.setEvents(events);

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" );
  output->cd();
  auto fr = doFit(make_pdf<EventList_type>(sig), events, eventsMC, MPS);

  auto ff = sig.fitFractions( fr->getErrorPropagator() );
  fr->addFractions(ff);
  fr->writeToFile( logFile );
  output->Close();
  return 0;
}
