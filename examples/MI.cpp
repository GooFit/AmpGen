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

#include "AmpGen/QcGenerator.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include <boost/algorithm/string/replace.hpp>
using namespace AmpGen;
using namespace std::complex_literals;


template <class PDF_TYPE, class PRIOR_TYPE> 
  void GenerateEvents( EventList& events
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& prior
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm )
{
  Generator<PRIOR_TYPE> signalGenerator( prior );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, events, nEvents );
}

template <class PDF_TYPE, class PRIOR_TYPE> 
  void GenerateCorrEvents( EventList& eventsSig
                       , EventList& eventsTag
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& priorSig
                       , PRIOR_TYPE& priorTag
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm 
		       )
{
  QcGenerator<PRIOR_TYPE> signalGenerator( priorSig, priorTag );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, eventsSig, eventsTag, nEvents );
}


Event makeEvent(double s01, double s02, double eps, EventType type){
TRandom3 rndm;
  rndm.SetSeed( 0 );
  gRandom = &rndm;
    PhaseSpace ps(type);

    Event evt =ps.makeEvent();
    double x = evt.s(0,1);
    double y = evt.s(0,2);
    double dx = x - s01;
    double dy = y - s02;
    double dr = std::pow(std::pow(dx, 2) + std::pow(dy, 2), 0.5);
  while ( dr > eps ){


    evt =ps.makeEvent();

    //INFO("s01, s02 = "<<evt.s(0, 1)<<" , "<<evt.s(0,2) << ", "<<s01<< " , "<<s02);
    x = evt.s(0,1);
    y = evt.s(0,2);
    dx = x - s01;
    dy = y - s02;
    dr = std::pow(std::pow(dx, 2) + std::pow(dy, 2), 0.5);

  }

    return evt;
    



}

EventList makeEvents(std::vector<double> s01, std::vector<double> s02, double eps, EventType type){
    EventList evts(type);
    std::vector<Event> evts_v;
    for (int i=0;i<s01.size();i++){
        Event evt = makeEvent(s01[i], s02[i], eps, type);
        //INFO("made evt "<<i<<" s01 = "<< evt.s(0,1));
        evts.push_back(evt);

    }
    INFO("Made "<<evts.size()<<" events");
    return evts;
}



double strongPhaseDiff(CoherentSum A, CoherentSum C, Event event){
    //auto evtT = event;
    //evtT.swap(1,2);
    auto z =  A.getValNoCache(event) * std::conj(C.getValNoCache(event));
    return M_PI - std::imag(std::log(z/std::abs(z))) ;
}

int getBin(Event event, size_t nBins, CoherentSum A, CoherentSum C){
    int b=0;
    double dd = strongPhaseDiff(A, C, event);
    double s01 = event.s(0, 1);
    double s02 = event.s(0, 2);
    double dd0 = M_PI - dd;
    double ddInv = M_PI + dd0;
    for (int i=1;i<nBins+1;i++){
//        double b1 = -M_PI + i * 2 * M_PI/nBins;
        //double b2 = b1 + 2 * M_PI/nBins;
        double b1 = (2 * M_PI/nBins) * (i - (3/2));
        double b2 = (2 * M_PI/nBins) * (i - (1/2));



        
        if (dd>b1 && dd<b2 && s01 > s02) b = i;
        if ((ddInv)>b1 && (ddInv)<b2 && s01 < s02) b = -(i);
    }
    return b;
}

std::map< int, EventList > binEvents(EventList events, size_t nBins, CoherentSum A, CoherentSum C){
    std::map<int, EventList > map;
    std::vector<int> bins;

    for (int j=0;j<nBins;j++){
        std::pair<int, EventList > pP(j+1, EventList(events.eventType()));
        std::pair<int, EventList > pM(-(j+1), EventList(events.eventType()));
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

complex_t cs_i(CoherentSum A, CoherentSum C, EventList events, int N){
    complex_t num = 0;
    double denA2 = 0;
    double denC2 = 0;
    for(int i=0;i<N;i++){
        Event evt = events[i]; 
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
    std::map<int, EventList > m = binEvents(events, nBins, A, C);
    for (int i=0;i<nBins;i++){
        EventList evtP = m[i+1];
        EventList evtM = m[-(i+1)];
        complex_t zP = cs_i(A, C, evtP, evtP.size());
        complex_t zM = cs_i(A, C, evtM, evtM.size());
        std::pair<int, complex_t> pP(i+1, zP);
        std::pair<int, complex_t> pM(-(i+1), zM);
        cs.insert(pP);
        cs.insert(pM);
    }
    return cs;
}




double fB(double f, double fbar, double x, double y, double c, double s){
    return f + fbar * (std::pow(x, 2) + std::pow(y, 2)) + 2 * std::pow(f*fbar, 0.5) *( c * x + s * y);
}

double fCP(double f, double fbar,  double c,int CP){
    return (f  + fbar - 2 *CP * std::pow(f * fbar, 0.5) * c)/2;
}

double fKspipi(double f1, double f2, double fbar1, double fbar2, double c1, double c2, double s1, double s2){
    return (f1 * fbar2 + fbar1 * f2 + 2 * std::pow(f1 * f2 * fbar1 * fbar2, 0.5) * (c1 * c2 - s1 * s2))/2;
}

double getF(int i, std::map<int, EventList> map){
    int sum = 0;
    for (auto p : map){
        sum += p.second.size();
    }
    return map[i].size()/sum;
}

std::map<int, double> getFs(std::map<int, EventList> map){
    std::map<int, double> m;
    for (auto p : map){
        double f = getF(p.first, map);
        std::pair<int, double> fPair(p.first, f);
        m.insert(fPair);
    }
    return m;
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

  int NIntBins = NamedParameter<int>("NIntBins", 1000);
  INFO("Using "<<NIntBins<<" integration events");
  EventType eventType(pNames);
  EventList mc =  Generator<>(eventType, &rndm).generate(NIntBins);
  CoherentSum A(eventType, MPS);
  bool conjHead = NamedParameter<bool>("conjHead", true);
  CoherentSum C(eventType.conj(conjHead), MPS);
  A.setEvents(mc);
  C.setEvents(mc);
  A.setMC(mc);
  C.setMC(mc);
  A.prepare();
  C.prepare();



    
    /*
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

*/
/*
    ProfileClock t_chi2;
    t_chi2.start();
    double chi2_1 = chi2();
    t_chi2.stop();
    INFO("t(chi2) = "<<t_chi2);
    MPS["D0{K*(892)+[BW]{K0S0,pi+},pi-}_Re"]->setCurrentFitVal(0);
    double chi2_2 = chi2();
    INFO("chi2_1 = "<<chi2_1);
    INFO("chi2_2 = "<<chi2_2);

    ProfileClock t_makeEvent;
    t_makeEvent.start();
    Event evt = makeEvent(1.5, 1.5, eps, eventType);
    t_makeEvent.stop();
    INFO("s01, s02 = "<<evt.s(0, 1)<<" , "<<evt.s(0,2)<<" took "<<t_makeEvent);
    f->Close();
    */
    double eps = NamedParameter<double>("eps", 0.1);
    std::string refBinningFile = NamedParameter<std::string>("refBinningFile", "ref_equal.root");
    bool remakeReferenceBins = NamedParameter<bool>("remakeRefBins", false);
    TFile *fRefEqual;
    if (remakeReferenceBins){
        fRefEqual = TFile::Open(refBinningFile.c_str(), "UPDATE");
    }
    else{


        fRefEqual = TFile::Open(refBinningFile.c_str());
    }
    INFO("Opened ref file to make true bins");
    fRefEqual->cd();
    std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
    std::vector<EventList> my_bins;
    for (int bin_idx=0;bin_idx<bins.size();bin_idx++){

        int bin = bins[bin_idx];
        INFO("bin = "<<bin);


        if (remakeReferenceBins){ 
            TTree * t_bin;
            
            if (bin>0) t_bin = (TTree*)fRefEqual->Get(("bin"+std::to_string(bin)).c_str() );
            if (bin<0) t_bin = (TTree*)fRefEqual->Get(("binm"+std::to_string(std::abs(bin))).c_str() );

            
            Float_t s01_bin_i, s02_bin_i;
            t_bin->SetBranchAddress("s01", &s01_bin_i);
            t_bin->SetBranchAddress("s02", &s02_bin_i);
            std::vector<double> s01_bin, s02_bin;
            for (int i=0;i<t_bin->GetEntries();i++){
                t_bin->GetEntry(i);
                s01_bin.push_back((double)s01_bin_i);
                s02_bin.push_back((double)s02_bin_i);
            }
            INFO("n(s01) = "<<s01_bin.size());
            EventList my_bin_bin = makeEvents(s01_bin, s02_bin, eps, eventType);
            my_bins.push_back(my_bin_bin);
            if (bin>0) my_bin_bin.tree(("myBin" + std::to_string(bin)).c_str())->Write(("myBin" + std::to_string(bin)).c_str());
            if (bin<0) my_bin_bin.tree(("myBinm" + std::to_string(std::abs(bin))).c_str())->Write(("myBinm" + std::to_string(std::abs(bin))).c_str());
        }
        else{
            EventList my_bin_bin;
    //       TTree * my_t_bin;
            if (bin>0) my_bin_bin = EventList(("ref_equal.root:myBin"+std::to_string(bin)).c_str(), eventType );
            if (bin<0) my_bin_bin = EventList(("ref_equal.root:myBinm"+std::to_string(std::abs(bin))).c_str(), eventType );
            my_bins.push_back(my_bin_bin);
        }
    }

    fRefEqual->Close();

   auto chi2_BESIII = [&A, &C, nBins, my_bins, &zBESIII, NIntBins](){
       double _chi2=0;
       (*(&A)).prepare();
       (*(&C)).prepare();
        int sign = NamedParameter<int>("cisi_sign", 1);

        std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
      // std::map<int, complex_t> cisi = cs(A, C, mc, nBins);


        #pragma omp parallel for reduction( +: _chi2)
        for (int bin_idx=0;bin_idx<bins.size();bin_idx++){
//           complex_t myZP = cisi[i+1];
///           complex_t myZM = cisi[-(i+1)];
            
           int actualBin = bins[bin_idx];
           EventList my_bin = my_bins[bin_idx];
           complex_t cisi = cs_i(A, C, my_bin, NIntBins);
            
           complex_t besiiiZ = (*(&zBESIII))[actualBin];
           _chi2 += std::pow( sign * std::real(cisi) - std::real(besiiiZ), 2) + std::pow(sign * std::imag(cisi) - std::imag(besiiiZ), 2);

       }
       return _chi2;
   };

    auto chi2_binning = [my_bins, &A, &C, nBins, NIntBins](){
        double _chi2 = 0;
        std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
        (*(&A)).prepare();
        (*(&C)).prepare();
        
        
        for (int bin_idx=0;bin_idx<bins.size();bin_idx++){


            int actualBin = bins[bin_idx];
            EventList my_bin = my_bins[bin_idx];

            int numEvents = NIntBins;
            if (my_bin.size()<numEvents) numEvents = my_bin.size();

            #pragma omp parallel for reduction( +: _chi2)
            for (int i=0;i<numEvents;i++){
                int myBin = getBin(my_bin[i], nBins, (*(&A)), (*(&C)));
                _chi2 += std::pow(myBin - actualBin, 2);
//int getBin(Event event, size_t nBins, CoherentSum A, CoherentSum C){
            }
        }
 
        return _chi2;


    };
//    auto mychi2 = chi2_binning;
    
    CoherentSum A_mag(eventType, MPS);
    EventList magEvents;

    PhaseSpace phsp(eventType,&rndm);
    size_t nEvents      = NamedParameter<size_t>     ("nMagEvents"  , 10000, "Total number of events to generate" );
    size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );

    TRandom3 rand;
    rand.SetSeed( seed + 934534 );

    GenerateEvents( magEvents, A, phsp , nEvents, blockSize, &rand );
    A_mag.setEvents(magEvents);

    int NInt = NamedParameter<int>("NInt", 1e7);
    EventList magMC =  Generator<>(eventType, &rndm).generate(NInt);

    A_mag.setMC(magMC);

    auto LLMagFit = make_likelihood(magEvents, A_mag);
    double chi2_weight = NamedParameter<double>("chi2_weight", 1);
    INFO("chi2_weight = "<<chi2_weight);
    auto fitMagAndBin = [chi2_binning, &LLMagFit, chi2_weight](){
        double combLL = chi2_weight *  chi2_binning() + (*(&LLMagFit)).getVal();
        return combLL;
    };
    auto fitMagAndBESIII = [chi2_BESIII, &LLMagFit, chi2_weight](){
        double combLL = chi2_weight *  chi2_BESIII() + (*(&LLMagFit)).getVal();
        return combLL;
    };



    double LLMagFit_PreFit = LLMagFit.getVal();
    double chi2_binning_PreFit = chi2_binning();
    double fitMagAndBin_PreFit = fitMagAndBin();
    INFO("LLMag = "<<LLMagFit_PreFit);
    INFO("chi2_binning = "<<chi2_binning_PreFit);
    INFO("fitMagAndBin = "<<fitMagAndBin_PreFit);




    if( NamedParameter<bool>("doFit", false)){
        if (NamedParameter<bool>("fitBinning", true)){
        Minimiser mini(fitMagAndBin, &MPS);
        auto res_fit = mini.doFit();
        double LLMagFit_PostFit = LLMagFit.getVal();
        double chi2_binning_PostFit = chi2_binning();
        double fitMagAndBin_PostFit = fitMagAndBin();

        double diff_LLMagFit = LLMagFit_PostFit - LLMagFit_PreFit;
        double diff_chi2_binning = chi2_binning_PostFit - chi2_binning_PreFit;
        double diff_fitMagAndBin = fitMagAndBin_PostFit - fitMagAndBin_PreFit;

        INFO("LLMag = "<<LLMagFit_PostFit);
        INFO("chi2_binning = "<<chi2_binning_PostFit);
        INFO("fitMagAndBin = "<<fitMagAndBin_PostFit);

        INFO("dLLMag = "<<diff_LLMagFit);
        INFO("dchi2_binning = "<<diff_chi2_binning);
        INFO("dfitMagAndBin = "<<diff_fitMagAndBin);

        FitResult * fr = new FitResult(mini);
        fr->writeToFile(logFile);
        }
        else if (NamedParameter<bool>("fitBESIII", false)){
            INFO("Fitting to given cisi from BESIII (2020)");

        Minimiser mini(fitMagAndBESIII, &MPS);
        mini.doFit();
        FitResult * fr = new FitResult(mini);
        fr->writeToFile(logFile);



        }
        else if (NamedParameter<bool>("fitMag", false)){
            INFO("Fitting magnitude only - gen + fit same model");
            Minimiser mini(LLMagFit, &MPS);
            mini.doFit();
            FitResult * fr = new FitResult(mini);
            fr->writeToFile(logFile);
        }
        
    }
   INFO("Writing bins + cisi");
  size_t NIntBinPlot = NamedParameter<size_t>("NIntBinPlot", 10000);
  EventList mcPlot =  Generator<>(eventType, &rndm).generate(NIntBinPlot);
  CoherentSum APlot(eventType, MPS);
  CoherentSum CPlot(eventType.conj(conjHead), MPS);
  APlot.setEvents(mcPlot);
  CPlot.setEvents(mcPlot);
  APlot.setMC(mcPlot);
  CPlot.setMC(mcPlot);
  APlot.prepare();
  CPlot.prepare();



    TFile* f = TFile::Open( output.c_str(), "RECREATE" );
    TNtuple * tup = new TNtuple( "Vals", "Vals", "aR:aI:s01:s02:dd:bin");
   TNtuple * tup_cs = new TNtuple("cisi", "cisi", "i:c:s");
    TNtuple * tup_BESIII = new TNtuple("cisi_BESIII", "cisi_BESIII", "i:c:s");
    INFO("Binning events for plotting");
   std::map<int, EventList > mPlot = binEvents(mcPlot, nBins, APlot, CPlot);
   int count=0;
   for (auto p : mPlot){
       int bin = p.first;
       EventList evts = p.second;
       for (auto evt:evts){
            double s01 = evt.s(0,1);
            double s02 = evt.s(0,2);
            double dd = strongPhaseDiff(APlot, CPlot, evt);
            complex_t a = APlot.getValNoCache(evt);
            Event evtT = evt;
            evtT.swap(1,2);
            complex_t aT = std::conj(APlot.getValNoCache(evtT));
            complex_t c = CPlot.getValNoCache(evt);
           // if (count*100/NIntBinPlot % 5 == 0) INFO("c = "<<c<<" aT = "<<aT<< " a = "<<a);
            count++;
            tup->Fill(std::real(a), std::imag(a), s01, s02, dd, bin);
       }
       complex_t z = cs_i( APlot, CPlot, evts, evts.size());

        tup_cs->Fill(bin, std::real(z), std::imag(z));
   }
   tup->Write("Vals");

    tup_cs->Write("cisi");

   INFO("writing BESIII cisi");


    for (int i=0;i<cBESIII.size();i++){
        int bin = i+1;
        tup_BESIII->Fill(bin, cBESIII[i], sBESIII[i]);
        tup_BESIII->Fill(-bin, cBESIII[i], -sBESIII[i]);
    }
    tup_BESIII->Write("cisi_BESIII");




/*
    std::map<int, complex_t> cisi = cs(A, C, mcPlot, nBins);
 
    for (auto p:cisi){
        int bin = p.first;
        complex_t z = p.second;
        tup_cs->Fill(bin, std::real(z), std::imag(z));
    }
    tup_cs->Write("cisi");
    */
    INFO("Closing file");

    f->Close();




   
    return 0;
}

/*
void GGSZ(){
    EventType KK("D0 K+ K-");
    EventType Ks("D0 K0S0 pi0");
    EventType Kppim("D0 K+ pi-");

    EventList sig_KK(eventType);
    EventList sig_Kspipi(eventType);
    EventList tag_KK(KK);
    EventList tag_Kppim(Kppim);
    EventList tag_Kspipi(Kspipi);


	
    PhaseSpace phspSig(eventType,&rndm);
    PhaseSpace phspTag_KK(KK,&rndm);
    pCorrelatedSum cs_KK(eventType, KK, MPS);
    pCorrelatedSum cs_Kppim(eventType, Kppim, MPS);
    pCorrelatedSum cs_Kspipi(eventType, eventType, MPS);
    INFO("Generating Events now!");

    size_t nCorrEvents      = NamedParameter<size_t>     ("nCorrEvents"  , 10000, "Total number of events to generate" );
    GenerateCorrEvents( sig_KK, tag_KK, cs_KK, phspSig_KK, phspTag_KK , nCorrEvents, blockSize, &rndm );
    GenerateCorrEvents( sig_Kppim, tag_Kppim, cs_Kppim, phspSig_Kppim, phspTag_Kppim , nCorrEvents, blockSize, &rndm );
    GenerateCorrEvents( sig_Kspipi, tag_Kspipi, cs_Kspipi, phspSig, phspTag , nCorrEvents, blockSize, &rndm );

    pCoherentSum psiBplus( eventType, MPS, 1, true, true);
    pCoherentSum psiBminus( eventType, MPS, -1, true, false);
    EventList BpEvents(eventType);
    EventList BmEvents(eventType);
    
    size_t nBEvents      = NamedParameter<size_t>     ("nBEvents"  , 10000, "Total number of events to generate" );
    GenerateEvents(BpEvents, psiBplus, phsp, nBEvents, blockSize, &rndm);
    GenerateEvents(BmEvents, psiBminus, phsp, nBEvents, blockSize, &rndm);

    int nBins = NamedParameter<int>("nBins", 8, "number of bins for MI");
    
    EventList mc =  Generator<>(eventType, &rndm).generate(NIntBins);
    CoherentSum A(eventType, MPS);
    bool conjHead = NamedParameter<bool>("conjHead", true);
    CoherentSum C(eventType.conj(conjHead), MPS);
    A.setEvents(mc);
    C.setEvents(mc);
    A.setMC(mc);
    C.setMC(mc);
    A.prepare();
    C.prepare();



    std::map<int, EventList> m_KK = binEvents(sig_KK, nBins, A, C);
    std::map<int, EventList> m_Kppim = binEvents(sig_Kppim, nBins, A, C);
    std::map<int, EventList> m_Kspipi_Sig = binEvents(sig_Kspipi, nBins, A, C);
    std::map<int, EventList> m_Kspipi_Tag = binEvents(tag_Kspipi, nBins, A, C);
    std::map<int, EventList> m_Bp = binEvents(BpEvents, nBins, A, C);
    std::map<int, EventList> m_Bm = binEvents(BmEvents, nBins, A, C);
    std::map<int, double> F_KK = getFs(m_KK);
    std::map<int, double> F_Kppim = getFs(m_Kppim);
    std::map<int, double> F_Bp = getFs(m_Bp);
    std::map<int, double> F_Bm = getFs(m_Bm);
    std::map<int, double> F_Kspipi_Sig = getFs(m_Kspipi_Sig);
    std::map<int, double> F_Kspipi_Tag = getFs(m_Kspipi_Tag);
    
    auto chi2_KK = [m_KK, MPS, F_KK, F_Kppim](){
        double _chi2 = 0;
        std
        for (auto p : m_KK){
            int bin = p.first;
            int binSign = bin/std::abs(bin);
            EventList evts = p.second;
            double c = MPS[("c" + std::to_string(std::abs(bin)).c_str() ]->mean();
            double s = binSign * MPS[("s" + std::to_string(std::abs(bin)).c_str() ]->mean();

            double Fi = F_Kppim[bin];
            double Fbari = F_Kppim[-bin];

            double fExp = fCP(Fi, Fbari, c, -1);
        




}


        }
    }
    







}
*/
