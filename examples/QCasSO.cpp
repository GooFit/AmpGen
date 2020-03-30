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
#include "AmpGen/ErrorPropagator.h"
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace AmpGen;

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data, EventList& mc, MinuitParameterSet& MPS, int nBins );
std::vector<std::string> makeBranches(EventType Type, std::string prefix);

std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps);
std::map<std::string, double > getPulls(std::map<std::string, std::vector<double> > fits, std::map<std::string, std::vector<double> > inits);
void writePulls(std::string fileName, std::map<std::string, double> pulls);

int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */

  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();

  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  

  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();

  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  size_t      nEvents  = NamedParameter<size_t>     ("nEvents"   , 10000       , "Number of events to fill in") ;
  int nBins = NamedParameter<int> ("nBins", 100, "number of bins for projection");

  bool QcGen2 = NamedParameter<bool>("QcGen2", false, "internal boolean - for new QcGenerator");
  bool NoQC = NamedParameter<bool>("noQC", false, "internal boolean - debug - no QC at all");

  bool DoFit = NamedParameter<bool>("doFit", true, "Do the fit");


   


  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

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
for( auto& tag : tags )
 {
    EventType signalType( pNames );
    auto tokens       = split(tag, ' ');
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType(); 
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");

    EventList sigEvents;
    EventList tagEvents;
    EventList sigMCEvents;
    EventList tagMCEvents;
    EventList UsigEvents;
    EventList UsigMCEvents;
        //EventList eventsMC = intFile == "" ? Generator<>(signalType, &rndm).generate(2e6) : EventList(intFile, signalType, GetGenPdf(true));
    EventList eventsMC;

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  

  std::map<std::string, std::vector<double> > inits = getParams(MPS);


    if (QcGen2){
        INFO("Loading QC Events");
      sigEvents = EventList(dataFile + ":Signal" , signalType);
      tagEvents = EventList(dataFile + ":Tag" , tagType);
//      eventsMC = EventList(intFile + ":Signal", signalType);

      eventsMC = EventList(intFile, signalType);
      

    }
    else if (NoQC){
      INFO("No QC");
      sigEvents = EventList(dataFile, signalType);
      eventsMC = EventList(intFile, signalType);
    }
    else{
    sigEvents = EventList(dataFile +":"+ tokens[0], signalType, Branches(sigBranches));

    tagEvents = EventList(dataFile +":"+ tokens[0], tagType, Branches(tagBranches));

      eventsMC = EventList(intFile + ":" +  tokens[0], signalType, Branches(sigBranches));

    }



  /* An EventType specifies the initial and final state particles as a vector that will be described by the fit. 
     It is typically loaded from the interface parameter EventType. */

  /* A CoherentSum is the typical amplitude to be used, that is some sum over quasi two-body contributions 
     weighted by an appropriate complex amplitude. The CoherentSum is generated from the couplings described 
     by a set of parameters (in a MinuitParameterSet), and an EventType, which matches these parameters 
     to a given final state and a set of data. A common set of rules can be matched to multiple final states, 
     i.e. to facilitate the analysis of coupled channels. 
     The CoherentSum is only appropriate for decays involving only (pseudo)scalars in the inital / final state, 
     otherwise the sum must also be over initial / final spin states. In this case, as PolarisedSum should be used. 
     See FitterWithPolarisation for an example of this use case.    
  */
  CoherentSum sig(signalType, MPS);
  
  /* Events are read in from ROOT files. If only the filename and the event type are specified, 
     the file is assumed to be in the specific format that is defined by the event type, 
     unless the branches to load are specified in the user options */

  
  /* Generate events to normalise the PDF with. This can also be loaded from a file, 
     which will be the case when efficiency variations are included. Default number of normalisation events 
     is 5 million. */
//  EventList eventsMC = intFile == "" ? Generator<>(signalType, &rndm).generate(2e6) : EventList(intFile, signalType, GetGenPdf(true));
  
  sig.setMC( eventsMC );

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  
  /* Do the fit and return the fit results, which can be written to the log and contains the 
     covariance matrix, fit parameters, and other observables such as fit fractions */
  FitResult* fr = doFit(make_likelihood(sigEvents, sig), sigEvents, eventsMC, MPS , nBins);
  /* Calculate the `fit fractions` using the signal model and the error propagator (i.e. 
     fit results + covariance matrix) of the fit result, and write them to a file. 
   */
  auto fitFractions = sig.fitFractions( fr->getErrorPropagator() ); 
  
  fr->addFractions( fitFractions );
  fr->writeToFile( logFile );
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(std::map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(logFile, pulls);
     

  output->cd();
  
  /* Write out the data plots. This also shows the first example of the named arguments 
     to functions, emulating python's behaviour in this area */

//  auto plots = sigEvents.makeDefaultProjections(Prefix("Data"), Bins(nBins));
 // for ( auto& plot : plots ) plot->Write();


    auto projections = signalType.defaultProjections(nBins);
    for( auto& projection : projections ){
      auto data_plot = projection(sigEvents);
      auto hist = projection.plot();
      for(unsigned i = 0 ; i != eventsMC.size(); ++i)
      {
        hist->Fill( projection( eventsMC[i] ), sig.prob( eventsMC[i]) );
      }
      hist->Scale( data_plot->Integral() / hist->Integral() );
      hist->SetName( (std::string("MC_")+hist->GetName()).c_str() );
      hist->Write();
      data_plot->Write();
      auto pull = (TH1D*) hist->Clone();
      pull->Add(data_plot, -1);
      pull->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
      pull->Write();
    }
    auto p2 = signalType.defaultProjections(nBins);
    for( unsigned i = 0 ; i != p2.size() -1; ++i )
    {
      for( unsigned j=i+1; j < p2.size(); ++j )
      {
        auto dalitz = Projection2D( p2[i], p2[j] );
        auto hdalitz = dalitz.plot();
        auto data_plot = sigEvents.makeProjection(dalitz);
        for( unsigned event = 0 ; event != eventsMC.size(); ++event )
        {
          auto pos = dalitz(eventsMC[event]);
          hdalitz->Fill( pos.first, pos.second, sig.prob( eventsMC[event]) );
        }
        hdalitz->Scale( data_plot->Integral() / hdalitz->Integral() );
        hdalitz->SetName( ( std::string("MC_") + hdalitz->GetName() ).c_str() );
        hdalitz->Write();

        data_plot->Write(); 
        auto pull2D = (TH2D*) hdalitz->Clone();
        pull2D->Add(data_plot, -1);

        pull2D->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
        pull2D->Write();
      }
    }


  
 
  output->Close();
}
}

template <typename likelihoodType>
FitResult* doFit( likelihoodType&& likelihood, EventList& data, EventList& mc, MinuitParameterSet& MPS, int nBins )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
  Minimiser mini( likelihood, &MPS );

  auto covar = mini.covMatrix();
  INFO("Printing Covariant matrix");
  covar.Print();
  auto covarFull = mini.covMatrixFull();
  //INFO("Printing full Covariant matrix");
  //covarFull.Print();
  mini.gradientTest();
  mini.doFit();

  FitResult* fr = new FitResult(mini);
  
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );

  /* Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
     The structure assumed the PDF is some SumPDF<eventListType, pdfType1, pdfType2,... >. */
  unsigned int counter = 1;
  for_each(likelihood.pdfs(), [&](auto& pdf){
    auto pfx = Prefix("Model_cat"+std::to_string(counter));
    auto mc_plot3 = mc.makeDefaultProjections(WeightFunction(pdf), Bins(nBins), pfx);
    for( auto& plot : mc_plot3 )
    {
      plot->Scale( ( data.integral() * pdf.getWeight() ) / plot->Integral() );
      plot->Write();
    }
    counter++;
  });

  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
  Chi2Estimator chi2( data, mc, likelihood, 15 );
  chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2.chi2(), chi2.nBins() );
  
  fr->print();
  return fr;
}
std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps){
  std::map<std::string, std::vector<double> > out;
  for (auto& param : mps){
    std::vector<double> vect = {param->mean(), param->err()};
    out[param->name()] = vect;
  }
  return out;
}




std::map<std::string, double > getPulls(std::map<std::string, std::vector<double> > fits, std::map<std::string, std::vector<double> > inits)  {
  std::string output = "";
 std::map<std::string, double> out;
 for(std::map<std::string, std::vector<double> >::iterator init = inits.begin(); init != inits.end(); ++init) {
  for(std::map<std::string, std::vector<double> >::iterator fit = fits.begin(); fit != fits.end(); ++fit) {
    std::string initName = init->first;
    std::vector<double> initParams = init->second;
    std::string fitName = fit->first;
    std::vector<double> fitParams = fit->second;
    if (initName == fitName){
        double pull = fitParams[0] - initParams[0];
        if (fitParams[1] != 0){
          pull /= fitParams[1];
        }
        out[fitName] = pull;
    }
  }
 }
  return out;

}


void writePulls(std::string fileName, std::map<std::string, double> pulls){
  std::ofstream outfile;
  outfile.open(fileName, std::ios_base::app);
      for(std::map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
        outfile<<"Pull "<<it->first<<" "<<it->second<<"\n";
       }

}



std::vector<std::string> makeBranches(EventType Type, std::string prefix){
  auto n = Type.finalStates().size();
  std::vector<std::string> branches;
  std::vector<std::string> varNames = {"PX", "PY", "PZ", "E"};
  for (long unsigned int i=0; i<n; i++){
    auto part = replaceAll(Type.finalStates()[i], "+", "p");
    part = replaceAll(part, "-", "m");
    for (auto varName : varNames){
      std::ostringstream stringStream;
      stringStream<<prefix<<part<<"_"<<varName;
      branches.push_back(stringStream.str());
    }
  }
  return branches;
}
