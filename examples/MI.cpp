#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/CombCorrLL.h"
#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/polyLASSO.h"
#include "AmpGen/ProfileClock.h"
#include <TMath.h>
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CombGamLL.h"
//#include <Math/IFunction.h>
#include <Math/Functor.h>
#include <TGraph.h>
#include <Minuit2/Minuit2Minimizer.h>
#include "TNtuple.h"
#include "AmpGen/PhaseCorrection.h"
#include <typeinfo>


#include <boost/algorithm/string/replace.hpp>
using namespace AmpGen;
using namespace std::complex_literals;


double strongPhaseDiff(CoherentSum A, CoherentSum C, Event event){
    //auto evtT = event;
    //evtT.swap(1,2);
    auto z =  A.getValNoCache(event) * std::conj(C.getValNoCache(event));
    return std::imag(std::log(z/std::abs(z)));
}

int getBin(Event event, size_t nBins, CoherentSum A, CoherentSum C){
    int b=0;
    double dd = strongPhaseDiff(A, C, event);
    double s01 = event.s(0, 1);
    double s02 = event.s(0, 2);

    for (int i=0;i<nBins;i++){
        double b1 = -M_PI + i * 2 * M_PI/nBins;
        double b2 = b1 + 2 * M_PI/nBins;
        if (dd>=b1 && dd<=b2 && s01 < s02) b = i+1;
        if ((-dd)>=b1 && (-dd)<=b2 && s01 > s02) b = -(i+1);
    }
    return b;
}

std::map< int, std::vector<Event> > binEvents(EventList events, size_t nBins, CoherentSum A, CoherentSum C){
    std::map<int, std::vector<Event> > map;
    std::vector<int> bins;

    for (int j=0;j<nBins;j++){
        std::pair<int, std::vector<Event> > pP(j+1, {});
        std::pair<int, std::vector<Event> > pM(-(j+1), {});
        map.insert(pP);
        map.insert(pM);


    }

    for (int j=0;j<events.size();j++){
        double dd = strongPhaseDiff(A, C, events[j]);
        double s01 = events[j].s(0, 1);
        double s02 = events[j].s(0, 2);
        int bin = getBin(events[j], nBins, A, C);

        map[bin].push_back(events[j]);
        //std::pair<Event, int> p(events[j], bin);
        //map.insert(p);
        //bins.push_back(bin);
    }
//    INFO("map size = "<<map.size());

   // INFO("bin size = "<<bins.size());
    return map;
}

complex_t cs_i(CoherentSum A, CoherentSum C, std::vector<Event> events){
    complex_t num = 0;
    double denA2 = 0;
    double denC2 = 0;
    for(auto  evt:events){ 
        num += A.getValNoCache(evt) * std::conj(C.getValNoCache(evt));
        denA2 += std::norm(A.getValNoCache(evt));
        denC2 += std::norm(C.getValNoCache(evt)); 
    }
    double den = std::pow(denA2 * denC2, 0.5);
    if (den!=0) return num/den;
    return 0;
}

std::map<int, complex_t> cs(CoherentSum A, CoherentSum C, EventList events, int nBins){
    std::map<int, complex_t> cs;
    A.transferParameters();
    C.transferParameters();
    std::map<int, std::vector<Event> > m = binEvents(events, nBins, A, C);
    for (int i=0;i<nBins;i++){
        std::vector<Event> evtP = m[i+1];
        std::vector<Event> evtM = m[-(i+1)];
        complex_t zP = cs_i(A, C, evtP);
        complex_t zM = cs_i(A, C, evtM);
        std::pair<int, complex_t> pP(i+1, zP);
        std::pair<int, complex_t> pM(-(i+1), zM);
        cs.insert(pP);
        cs.insert(pM);
    }
    return cs;
}




int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );
  //OptionsParser::setArgs( argc, argv, "Toy simulation for Quantum Correlated Ψ(3770) decays");
  /* */
  //auto time_wall = std::chrono::high_resolution_clock::now();
  //auto time      = std::clock();
  size_t hwt = std::thread::hardware_concurrency();
  size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
  //double luminosity   = NamedParameter<double>("Luminosity"  , 818.3       , "Luminosity to generate. Defaults to CLEO-c integrated luminosity.");
  //size_t nEvents      = NamedParameter<size_t>("nEvents"     , 0           , "Can also generate a fixed number of events per tag, if unspecified use the CLEO-c integrated luminosity.");
  size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
  //bool   poissonYield = NamedParameter<bool  >("PoissonYield", true        , "Flag to include Poisson fluctuations in expected yields (only if nEvents is not specified)");
  //bool   noQCFit      = NamedParameter<bool  >("noQCFit"     , false       , "Treat Signal and Tag as uncorrelated and fit the data individually");
  double crossSection = NamedParameter<double>("CrossSection", 3.26 * 1000 , "Cross section for e⁺e⁻ → Ψ(3770) → DD'");
  std::string output  = NamedParameter<std::string>("Output" , "ToyMC.root", "File containing output events"); 
  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();  

  auto BTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();

  bool m_debug        = NamedParameter<bool>("Debug", false, "Debug QcFitter output");
    bool doDebugNorm  = NamedParameter<bool>("doDebugNorm", false, "Debug the normalisation of the pdf");
    int nBins = NamedParameter<int>("nBins", 8, "number of bins for MI");
    int fBins = NamedParameter<int>("fBins", 1, "fraction of nEvents to be used as bins for projections");
    int nFits = NamedParameter<int>("nFits", 4, "number of repeats of mini.doFits() for debug purposes!");
    bool doProjections = NamedParameter<bool>("doProjections", true);
    bool doScan = NamedParameter<bool>("doScan", true);
    bool doPCorrSum = NamedParameter<bool>("doPCorrSum", false);
//    bool doCombFit = NamedParameter<bool>("doCombFit", false, "Do a combined fit of 3 tags - at the moment this is hard coded for now");

    auto NIntMods = NamedParameter<std::string >("NIntMods", std::string("1"), "Number of modifiers for NInt - go in ascending order").getVector();
  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
  bool doCombFit      = NamedParameter<bool>("doCombFit", true, "Do combined fit");
  bool doTagFit      = NamedParameter<bool>("doTagFit", true, "Do fit for each tag");
  int  maxAttempts = NamedParameter<int>("maxAttempts", 5, "Max attempts to get a valid minimum from Minuit2");
  bool QcGen2 = NamedParameter<bool>("QcGen2", false, "internal boolean - for new QcGenerator");
  bool doFit = NamedParameter<bool>("doFit", true, "Do the fit");
  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);
  if (intFile == ""){

  }
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  std::vector<std::string> varNames = {"E", "PX", "PY", "PZ"};
  //auto yc = DTYieldCalculator(crossSection);
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
//    add_CP_conjugate(MPS);
      AddCPConjugate(MPS);
  }

    std::vector<double> cBESIII = {0.708,0.671,0.001,-0.602, -0.965, -0.554, 0.046, 0.403}; 
    std::vector<double> sBESIII = {0.128, 0.341, 0.893, 0.723, 0.020, -0.589, -0.686, -0.474};
    std::map<int, complex_t> zBESIII = {};
    for (int i=0;i<cBESIII.size();i++){
        complex_t zP(cBESIII[i], sBESIII[i]);
        complex_t zM(cBESIII[i], -sBESIII[i]);
        std::pair<int, complex_t> pP(i+1, zP);
        std::pair<int, complex_t> pM(-(i+1), zM);
        zBESIII.insert(pP);
        zBESIII.insert(pM);
    }

  int NInt = NamedParameter<int>("NInt", 1e7);
  INFO("Using "<<NInt<<" integration events");
  EventType eventType(pNames);
  EventList mc =  Generator<>(eventType, &rndm).generate(NInt);
  CoherentSum A(eventType, MPS);
  CoherentSum C(eventType.conj(true), MPS);
  A.setEvents(mc);
  C.setEvents(mc);
  A.setMC(mc);
  C.setMC(mc);
  A.prepare();
  C.prepare();



    TFile* f = TFile::Open( output.c_str(), "RECREATE" );
    mc.tree("Events")->Write("Events");
   TNtuple * tup = new TNtuple( "Vals", "Vals", "aR:aI:s01:s02:dd:bin");
   std::map<int, std::vector<Event> > m = binEvents(mc, nBins, A, C);
   for (auto p : m){
       int bin = p.first;
       std::vector<Event> evts = p.second;
       for (auto evt:evts){
            double s01 = evt.s(0,1);
            double s02 = evt.s(0,2);
            double dd = strongPhaseDiff(A, C, evt);
            complex_t a = A.getValNoCache(evt);
            Event evtT = evt;
            evtT.swap(1,2);
            complex_t c = std::conj(A.getValNoCache(evtT));
            tup->Fill(std::real(a), std::imag(a), s01, s02, dd, bin);
       }
   }
   tup->Write("Vals");

    std::map<int, complex_t> cisi = cs(A, C, mc, nBins);
    TNtuple * tup_cs = new TNtuple("cisi", "cisi", "i:c:s");
    for (auto p:cisi){
        int bin = p.first;
        complex_t z = p.second;
        tup_cs->Fill(bin, std::real(z), std::imag(z));
    }
    tup_cs->Write("cisi");
    TNtuple * tup_BESIII = new TNtuple("cisi_BESIII", "cisi_BESIII", "i:c:s");
    for (int i=0;i<cBESIII.size();i++){
        int bin = i+1;
        tup_BESIII->Fill(bin, cBESIII[i], sBESIII[i]);
    }
    tup_BESIII->Write("cisi_BESIII");

   f->Close();

   auto chi2 = [&A, &C, nBins, mc, &zBESIII](){
       double _chi2=0;
       (*(&A)).transferParameters();
       (*(&C)).transferParameters();


       std::map<int, complex_t> cisi = cs(A, C, mc, nBins);

       for (int i=0;i<nBins;i++){
           complex_t myZP = cisi[i+1];
           complex_t myZM = cisi[-(i+1)];
           complex_t besiiiZP = (*(&zBESIII))[i+1];
           complex_t besiiiZM = (*(&zBESIII))[-(i+1)];
           double _chi2P = std::pow(std::real(myZP) - std::real(besiiiZP), 2) +std::pow(std::imag(myZP) - std::imag(besiiiZP), 2) ;
           double _chi2M = std::pow(std::real(myZM) - std::real(besiiiZM), 2) +std::pow(std::imag(myZM) - std::imag(besiiiZM), 2) ;

           _chi2 += _chi2P + _chi2M;

       }
       return _chi2;
   };

    ProfileClock t_chi2;
    t_chi2.start();
    double chi2_1 = chi2();
    t_chi2.stop();
    INFO("t(chi2) = "<<t_chi2);
    MPS["D0{K*(892)+[BW]{K0S0,pi+},pi-}_Re"]->setCurrentFitVal(0);
    double chi2_2 = chi2();
    INFO("chi2_1 = "<<chi2_1);
    INFO("chi2_2 = "<<chi2_2);



   Minimiser mini(chi2, &MPS);
   auto res_fit = mini.doFit();
   //FitResult * fr = new FitResult(mini);
   //fr->writeToFile(logFile);
    return 0;
}