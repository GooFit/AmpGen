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
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/AddCPConjugate.h"
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

  const auto datasets        = NamedParameter<std::string>("Datasets","",
      "List of data/simulated samples to fit, in the format \
      \033[3m data[0] sim[0] data[1] sim[1] ... \033[0m. \nIf a simulated sample is specified FLAT, uniformly generated phase-space events are used for integrals ").getVector();

  const std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log",
      "Name of the output log file");

  const std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root",
      "Name of the output plot file");
  
  const bool add_conj = NamedParameter<bool>("AddConj", false );

  std::vector<EventList_type>             data;
  std::vector<EventList_type>              mcs;

  #ifdef _OPENMP
    size_t hwThreads = std::thread::hardware_concurrency();
    size_t usThreads = NamedParameter<size_t>( "nCores", hwThreads, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << usThreads  << " / " << hwThreads << " threads" );
    omp_set_num_threads(usThreads);
    omp_set_dynamic(0);
  #endif

  INFO("Output : " << logFile << " plots = " << plotFile );

  for(size_t i=0;i < datasets.size() ; i+=2 ){
    data.emplace_back( datasets[i] );
    if( datasets[i+1] == "FLAT" ) mcs.emplace_back(Generator<>(data.rbegin()->eventType() ).generate(1e6));
    else mcs.emplace_back( datasets[i+1], data.rbegin()->eventType(), GetGenPdf(true) );
  }

  std::vector<PolarisedSum>           fcs(data.size());
  std::vector<SumPDF<EventList_type, PolarisedSum&>> pdfs;

  pdfs.reserve(data.size());

  SimFit totalLL;
  MinuitParameterSet mps;
  mps.loadFromStream();
  if( add_conj ) AddCPConjugate(mps);
  for(size_t i = 0; i < data.size(); ++i){
    fcs[i] = PolarisedSum(data[i].eventType(), mps);
    pdfs.emplace_back( make_pdf<EventList_type>(fcs[i]) );
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
    for( auto proj : data[i].eventType().defaultProjections() )
    {
      proj(data[i], PlotOptions::Prefix("Data"+std::to_string(i)), PlotOptions::AutoWrite() );
      proj(mcs[i]  , pdfs[i].componentEvaluator(&mcs[i]), PlotOptions::Prefix("pdf"+std::to_string(i)), PlotOptions::Norm(data.size()), PlotOptions::AutoWrite() );
    }    
  }
  output_plots->Close();
}
