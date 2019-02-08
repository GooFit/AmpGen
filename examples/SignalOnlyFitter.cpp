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
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"

#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace AmpGen;

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data, EventList& mc, MinuitParameterSet& MPS );

int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  const std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  const std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", "Name of the output log file");
  const std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  auto evtType_particles     = NamedParameter<std::string>("EventType" , ""          , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  
  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( evtType_particles.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  [[maybe_unused]]
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  const size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
  
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  /* A MinuitParameterSet is (unsurprisingly) a set of fit parameters, and can be loaded from 
     the parsed options. For historical reasons, this is referred to as loading it from a "Stream" */
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  /* An EventType specifies the initial and final state particles as a vector that will be described by the fit. 
     It is typically loaded from the interface parameter EventType. */
  EventType evtType(evtType_particles);
  
  /* A CoherentSum is the typical amplitude to be used, that is some sum over quasi two-body contributions 
     weighted by an appropriate complex amplitude. The CoherentSum is generated from the couplings described 
     by a set of parameters (in a MinuitParameterSet), and an EventType, which matches these parameters 
     to a given final state and a set of data. A common set of rules can be matched to multiple final states, 
     i.e. to facilitate the analysis of coupled channels. */
  CoherentSum sig(evtType, MPS);
  
  /* Events are read in from ROOT files. If only the filename and the event type are specified, 
     the file is assumed to be in the specific format that is defined by the event type. */
  EventList events(dataFile, evtType);
  
  /* Generate events to normalise the PDF with. This can also be loaded from a file, 
     which will be the case when efficiency variations are included. Default number of normalisation events 
     is 5 million. */
  EventList eventsMC = Generator<>(evtType, &rndm).generate(5e6);
  
  sig.setMC( eventsMC );

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  
  /* Do the fit and return the fit results, which can be written to the log and contains the 
     covariance matrix, fit parameters, and other observables such as fit fractions */
  FitResult* fr = doFit(make_pdf(sig), events, eventsMC, MPS );
  fr->writeToFile( logFile );
  output->cd();
  
  /* Write out the data plots. This also shows the first example of the named arguments 
     to functions, emulating python's behaviour in this area */

  auto plots = events.makeDefaultPlots(Prefix("Data"), Bins(100));
  for ( auto& plot : plots ) plot->Write();

  output->Close();
}

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data, EventList& mc, MinuitParameterSet& MPS )
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

  /* Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
     The structure assumed the PDF is some SumPDF<T1,T2,...>. */
  unsigned int counter = 1;
  for_each(pdf.m_pdfs, [&]( auto& f ){
    std::function<double(const Event&)> FCN_sig = [&](const Event& evt){ return f.prob_unnormalised(evt) ; };
    auto mc_plot3 = mc.makePlots(mc.eventType().defaultProjections(100)
                                , WeightFunction(f)
                                , Prefix("Model_cat"+std::to_string(counter)));
    for( auto& plot : mc_plot3 )
    {
      plot->Scale( ( data.integral() * f.getWeight() ) / plot->Integral() );
      plot->Write();
    }
    counter++;
  } );

  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
  Chi2Estimator chi2( data, mc, pdf, 15 );
  fr->addChi2( chi2.chi2(), chi2.nBins() );

  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );
  fr->print();
  return fr;
}
