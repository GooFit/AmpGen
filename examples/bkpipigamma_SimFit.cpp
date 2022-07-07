#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/Projection.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

using namespace AmpGen;

int main(int argc , char* argv[] ){
  OptionsParser::setArgs( argc, argv );
  
  const auto pEventType      = NamedParameter<std::string>("EventType", std::vector<std::string>(), 
        "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector();
  
  const auto datasets        = NamedParameter<std::string>("Datasets","", 
      "List of data/simulated samples to fit, in the format \
      \033[3m data[0] sim[0] data[1] sim[1] ... \033[0m. \nIf a simulated sample is specified FLAT, uniformly generated phase-space events are used for integrals ").getVector();
  
  const std::string weightbr = NamedParameter<std::string>("Weight"    , ""        , "Name of the weights branch (sweights)." );

  auto databNames = NamedParameter<std::string>("Data_Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();

  auto simbNames = NamedParameter<std::string>("Sim_Branches", std::vector<std::string>()
              ,"List of branch names in the integration sample, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();


  const std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", 
      "Name of the output log file");
  
  const std::string plotFile = NamedParameter<std::string>("Plots"     , "", 
      "Name of the output plot file");

  const size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );

  const bool  perturb  = NamedParameter<bool> ("Perturb"   , 1           , "Flag to randomise starting parameters.");

  
  #ifdef _OPENMP 
    size_t hwThreads = std::thread::hardware_concurrency();
    size_t usThreads = NamedParameter<size_t>( "nCores", hwThreads, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << usThreads  << " / " << hwThreads << " threads" );
    omp_set_num_threads(usThreads);
    omp_set_dynamic(0);
  #endif
  
  INFO("Output : " << logFile << " plots = " << plotFile );

  // Check that we have an event type
  if( pEventType.size() == 0 ) FATAL("Must specify event format as EventType \033[3m parent daughter1 daughter2 ... \033[0m in options");
  
  // Set the event type
  const EventType eventType  = EventType( pEventType );

  // Load the data and simulated data event lists:
  std::vector<EventList>             data;
  std::vector<EventList>              mcs; 

  for(size_t i=0;i < datasets.size() ; i+=2 ){
    if (weightbr ==  ""){
      // Load data sets without reading the weight branch;  (real data == no GenPdf to be read)
      data.emplace_back( datasets[i], eventType, GetGenPdf(false), Branches(databNames));
    }
    else {
      // Load data sets reading the weight branch;  (real data == no GenPdf to be read)
      data.emplace_back( datasets[i], eventType, GetGenPdf(false), Branches(databNames), WeightBranch(weightbr));
    }
    if( datasets[i+1] == "FLAT" ){
      // Generate data sets for normalisation (NO ACCEPTANCE)
      mcs.emplace_back(Generator<>(eventType).generate(5e6));
      INFO("Generated: 5e6 events for integrals" );
    }
    else {
      // Load simulated data sets for normalisation:
      mcs.emplace_back( datasets[i+1], eventType, GetGenPdf(true),Branches(simbNames) );
    }
  }
  INFO( "Data events: " << data.size() );
  INFO( "MC events  : " << mcs.size() );


  // Build the vector of PDFs
  std::vector<PolarisedSum>  fcs(data.size()); 
  std::vector<SumPDF<EventList, PolarisedSum&>> pdfs; 
  
  pdfs.reserve(data.size()); 
  
  // Make the total likelihood and read starting values for the parameters
  SimFit totalLL; 
  MinuitParameterSet mps;
  mps.loadFromStream();

  // randomise the initial values (always done for blind params, only if perturb is ON for the rest)
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;
  for ( auto& param : mps ) {
    // for lambda gamma:  uniformly distributed in [0,1]
    if (param-> isBlind()) {
      param->setCurrentFitVal( rndm.Rndm() );
      continue;
    }
    // for the rest, perturb with a gaussian, (mean+/-  3*error)			     
    if ( perturb ) {
      if ( !param->isFree() ) continue;
      param->setCurrentFitVal( rndm.Gaus( param->mean(), 3.*param->err() ) );
    }
  }

  // Load data in Sim PDF
  for(size_t i = 0; i < data.size(); ++i){
    fcs[i] = PolarisedSum(eventType, mps);
    pdfs.emplace_back( make_pdf(fcs[i]) );
    pdfs[i].setEvents(data[i]);
    auto& mc = mcs[i];
    for_each( pdfs[i].pdfs(), [&mc](auto& pdf){pdf.setMC(mc);});  
    totalLL.add( pdfs[i] );
  }

  // Declare the minimiser:
  Minimiser mini( totalLL, &mps );
  mini.doFit();
  FitResult* fr = new FitResult(mini);
  /*  
  // Estimate the chi2 
  Chi2Estimator* chi2;
  float GChi2  = 0.;
  float GnBins = 0.;
  for( size_t i = 0 ; i < data.size(); ++i ){
    INFO("Computing Chi2 for sample: " << i << " ...");
    chi2 = new Chi2Estimator( data[i], mcs[i], fcs[i], 5 );
    GChi2 += chi2->chi2();
    GnBins += chi2->nBins();
  }
  fr->addChi2( GChi2, GnBins );
  */
  // Calculate the fit fractions
  auto fitFractions = fcs[0].fitFractions( fr->getErrorPropagator() ); 
  fr->addFractions( fitFractions );
  fr->writeToFile( logFile );
  fr->print();


  //ATTEMPT AT: Save MC events with per PDF component weight
  //int i = 0;
  //TFile* output_plots = TFile::Open( plotFile.c_str(), "RECREATE");
  //perAmplitudePlot(mcs[i], fcs[i]);

/*
  // Make Default projections
  if (plotFile != ""){
    TFile* output_plots = TFile::Open( plotFile.c_str(), "RECREATE");
    for( size_t i = 0 ; i < data.size(); ++i )
      {
	INFO("Making figures for sample: " << i << " ...");
	auto dataPlots = data[i].makeDefaultProjections( Prefix("Data_"+std::to_string(i))); 
	for( auto& p : dataPlots ) p->Write();
	size_t counter = 0;
	for_each(pdfs[i].pdfs(), [&]( auto& f ){
	    auto mc_plots = mcs[i].makeDefaultProjections(WeightFunction(f), 
							  Prefix("Model_sample_"+std::to_string(i)+"_cat"+std::to_string(counter)));
	    for( auto& plot : mc_plots )
	      {
		plot->Scale( ( data[i].integral() * f.getWeight() ) / plot->Integral() );
		plot->Write();
	      }
	    counter++;
	  } );
      }
    output_plots->Close();
    INFO("Plotting done.");
  }
*/
}
