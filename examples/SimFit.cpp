#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include <omp.h>

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
  
  const std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", 
      "Name of the output log file");
  
  const std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", 
      "Name of the output plot file");

  std::vector<EventList>             data;
  std::vector<EventList>              mcs; 
  
  #ifdef _OPENMP 
    size_t hwThreads = std::thread::hardware_concurrency();
    size_t usThreads = NamedParameter<size_t>( "nCores", hwThreads, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << usThreads  << " / " << hwThreads << " threads" );
    omp_set_num_threads(usThreads);
    omp_set_dynamic(0);
  #endif
  
  INFO("Output : " << logFile << " plots = " << plotFile );
  
  if( pEventType.size() == 0 ) FATAL("Must specify event format as EventType \033[3m parent daughter1 daughter2 ... \033[0m in options");

  const EventType eventType  = EventType( pEventType );

  for(size_t i=0;i < datasets.size() ; i+=2 ){
    data.emplace_back( datasets[i], eventType );
    if( datasets[i+1] == "FLAT" ) mcs.emplace_back(Generator<>(eventType).generate(1e6));
    else mcs.emplace_back( datasets[i+1], eventType, GetGenPdf(true) );
  }

  std::vector<PolarisedSum>           fcs(data.size()); 
  std::vector<SumPDF<EventList, PolarisedSum&>> pdfs; 
  
  pdfs.reserve(data.size()); 

  SimFit totalLL; 
  MinuitParameterSet mps;
  mps.loadFromStream();
  for(size_t i = 0; i < data.size(); ++i){
    fcs[i] = PolarisedSum(eventType, mps);
    pdfs.emplace_back( make_pdf(fcs[i]) );
    pdfs[i].setEvents(data[i]);
    auto& mc = mcs[i];
    for_each( pdfs[i].pdfs(), [&mc](auto& pdf){pdf.setMC(mc);});  
    totalLL.add( pdfs[i] );
  } 
  Minimiser mini( totalLL, &mps );
  mini.doFit();  
  FitResult(mini).writeToFile("Fitter.log");
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
}
