#include <TH1.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Factory.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"

#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#if ENABLE_AVX2
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include "TFile.h"
#include "TRandom3.h"

using namespace AmpGen;

std::vector<ThreeBodyCalculator> threeBodyCalculators( MinuitParameterSet& mps )
{
  std::vector<std::string> threeBodiesToIntegrate = NamedParameter<std::string>( "ThreeBodiesToIntegrate" ).getVector();
  std::vector<ThreeBodyCalculator> calculators;
  for ( auto& v : threeBodiesToIntegrate ) calculators.emplace_back( v, mps );
  return calculators;
}

void randomizeStartingPoint( MinuitParameterSet& MPS, TRandom3& rand, bool SplineOnly = false )
{
  double range = 5;
  for (auto& param : MPS) {
    if ( ! param->isFree() == 0 ) continue;
    if ( SplineOnly && param->name().find( "::Spline::" ) == std::string::npos ) continue;
    range = param->maxInit() - param->minInit();
    param->setInit( range * rand.Rndm() + param->meanInit() );
    std::cout << *param << std::endl;
  }
}

template <typename SIGPDF>
void addExtendedTerms( Minimiser& mini, SIGPDF& pdf, MinuitParameterSet& mps )
{
  std::vector<std::string> llConfigs = NamedParameter<std::string>( "LLExtend" ).getVector();

  for ( auto& ll_config : llConfigs ) {
    auto ll_name = split( ll_config, ' ' )[0];
    auto ll_term = Factory<IExtendLikelihood>::get( ll_name );
    if ( ll_term != nullptr ) {
      ll_term->configure( ll_config, pdf, mps );
      mini.addExtendedTerm( ll_term );
    } else {
      ERROR( "LL term : " << ll_name << " not recognised" );
    }
  }
}

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data, EventList& mc, MinuitParameterSet& MPS )
{
  INFO( "Type = " << type_string<PDF>() );
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  pdf.setEvents( data );

  Minimiser mini( pdf, &MPS );
  //addExtendedTerms( mini, std::get<0>( pdf.pdfs() ), MPS );
  auto threeBodyShapes     = threeBodyCalculators( MPS );
  unsigned int updateWidth = NamedParameter<unsigned int>( "UpdateWidth", 0 );
  if ( updateWidth ) {
    for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
  }
  unsigned int nIterations            = NamedParameter<unsigned int>( "nIterations", 0 );
  std::vector<std::string> SlowParams = NamedParameter<std::string>( "Release", "" ).getVector();
  std::vector<MinuitParameter*> slowParamPtrs;
  if ( nIterations != 0 ) {
    for ( auto& param : SlowParams ) {
      auto it      = MPS.find( param );
      if ( it != nullptr ) {
        slowParamPtrs.push_back( it );
        it->fix();
      } else {
        WARNING( "Trying to release non-existent parameter: " << param );
      }
    }
  }
  INFO( "Fitting PDF with " << pdf.nPDFs() << " components, iterating " << " " << nIterations + 1 << " times" );
  for ( unsigned int iteration = 0; iteration < nIterations + 1; ++iteration ) {
    mini.doFit();
    if ( iteration == 0 && nIterations != 0 ) {
      for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
      for ( auto& param : slowParamPtrs ) param->setFree(); /// release the parameter ///
    }
  }

  FitResult* fr          = new FitResult( mini );
  bool makePlots = NamedParameter<bool>( "MakePlots", true );

  if ( makePlots ) {
    auto ep = fr->getErrorPropagator();

    unsigned int counter = 1;
    for_each( pdf.pdfs(), [&]( auto& f ) {
        auto tStartIntegral2 = std::chrono::high_resolution_clock::now();
        auto mc_plot3 = mc.makeProjections( mc.eventType().defaultProjections(100), WeightFunction(f), PlotOptions::Prefix("tMC_Category"+std::to_string(counter) ) );
        auto tEndIntegral2   = std::chrono::high_resolution_clock::now();
        double t2            = std::chrono::duration<double, std::milli>( tEndIntegral2 - tStartIntegral2 ).count();
        INFO( "Time for plots = " << t2 );

        for ( auto& plot : mc_plot3 ) {
        plot->Scale( ( data.integral() * f.getWeight() ) / plot->Integral() );
        plot->Write();
        }
        counter++;
        } );
  }
  Chi2Estimator chi2( data, mc, pdf, MinEvents(15) );
  fr->addChi2( chi2.chi2(), chi2.nBins() );

  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );
  fr->print();
  return fr;
}

int main( int argc, char* argv[] )
{
  OptionsParser::setArgs( argc, argv );

  const std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  const std::string mcFile   = NamedParameter<std::string>("SimSample" , ""          , "Name of file containing normalisation sample." );
  const std::string fmcFile  = NamedParameter<std::string>("FlatMC"    , ""          , "Name of file containing events for computing physics integrals" );
  const std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", "Name of the output log file");
  const std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");

  [[maybe_unused]]
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );

  const size_t      NBins    = NamedParameter<size_t>     ("nBins"     , 100         , "Number of bins used for plotting.");
  const bool        perturb  = NamedParameter<bool>       ("Perturb"   , 0           , "Flag to randomise starting parameters.");
  const size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );

  std::vector<std::string> evtType_particles =  NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  INFO( "Output : " << logFile << " plots = " << plotFile );
  
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if ( dataFile == "" ) {
    ERROR( "No input data selected" );
    return -1;
  }
  if ( mcFile == "" ) WARNING( "No input simulation selected; using PHSP sample" );

  if ( perturb ) {
    for ( auto& param : MPS ) {
      if ( !param->isFree() ) continue;
      param->setCurrentFitVal( rndm.Gaus( param->mean(), param->err() ) );
    }
  }
  EventType evtType( evtType_particles );
  INFO( "Signal  = " << evtType << " OS = " << evtType.conj(true) );
  CoherentSum   sig( evtType, MPS );
  IncoherentSum bkg( evtType, MPS, "Inco" );
  CoherentSum misID( evtType.conj(true) , MPS );

  INFO( "fPDF = " << MPS["fPDF"]->mean() << " "
      << "fComb = " << MPS["fComb"]->mean() << " "
      << "fMisID = " << MPS["fMisID"]->mean() );

  sig.setWeight( MPS["fPDF"] );
  bkg.setWeight( MPS["fComb"] );
  misID.setWeight( MPS["fMisID"] );

  const std::string cut         = NamedParameter<std::string>( "Cut", "1" );
  const std::string simCut      = NamedParameter<std::string>( "SimCut", "1" );
  bool BAR                      = NamedParameter<bool>("Bar",false);

  EventList events( dataFile, !BAR ? evtType : evtType.conj() , Filter(cut) );
  EventList eventsMC = mcFile == "" ? EventList( evtType) : EventList( mcFile, !BAR ? evtType : evtType.conj(), Filter(simCut) ) ;
  
    auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
  if( NamedParameter<std::string>("Units", "GeV").getVal()  == "MeV") {
    INFO("Changing units from MeV -> GeV");
    events.transform( scale_transform );
  }
    eventsMC.transform( scale_transform );
  
  INFO( "Data events: " << events.size() );  
  INFO( "MC events  : " << eventsMC.size() );
  if( mcFile == "" ){
    Generator<>( evtType, &rndm ).fillEventListPhaseSpace( eventsMC, 5e6 );
    INFO("Generated: " << eventsMC.size() << " events for integrals" );
  }
  
  sig.setMC( eventsMC );
  bkg.setMC( eventsMC );
  misID.setMC( eventsMC );

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" );
  output->cd();
  FitResult* fr = nullptr;

  if ( MPS["fPDF"]->mean() == 1.0 )
    fr = doFit( make_pdf( sig ), events, eventsMC, MPS );

  else if ( MPS["fPDF"]->mean() == 0.0 && MPS["fComb"]->mean() == 1.0 )
    fr = doFit( make_pdf( bkg ), events, eventsMC, MPS );

  else if ( MPS["fMisID"]->mean() == 0.0 )
    fr = doFit( make_pdf( sig, bkg ), events, eventsMC, MPS );

  else if ( MPS["fMisID"]->mean() != 0.0 )
    fr = doFit( make_pdf( sig, bkg, misID ), events, eventsMC, MPS );

  if ( fr == nullptr ) {
    ERROR( "Fit fails" );
    return -1;
  }

  INFO( "Completed fit; calculating additional observables" );

//  EventList* fmc               = &eventsMC;
//  if ( fmcFile != mcFile ) fmc = new EventList( fmcFile, evtType );
//  EventList& flatMC            = *fmc;
//  sig.reset( true ); //// reset PDFs to ensure correct cache state
//  sig.setMC( flatMC );
//  sig.prepare();
  if ( MPS["fComb"]->mean() != 1 )
    fr->addFractions( sig.fitFractions( fr->getErrorPropagator() ) );
//  else
//    fr->addFractions( bkg.fitFractions( fr->getErrorPropagator() ) );

  fr->writeToFile( logFile );
  output->cd();
  auto plots = events.makeDefaultProjections( PlotOptions::Prefix( "Data_" ), PlotOptions::Bins( NBins ) );
  for ( auto& plot : plots ) plot->Write();

  output->Write();
  output->Close();

  INFO( "Finalising output" );
  return 0;
}
