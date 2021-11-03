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
#include "AmpGen/ProgressBar.h"
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

double minS(Event event, EventList events, double fracNevents=1){
    double dS = 100;
    int nEvents = events.size() * fracNevents;
//    INFO("Going through "<<nEvents);
    if (nEvents == 0) nEvents = 1;

    std::vector<double> dSs(nEvents);
    #pragma omp parallel for
    for (int i=0;i<nEvents;i++){
        Event evt = events[i];
        double dx = event.s(0,1) - evt.s(0,1);
        double dy = event.s(0,2) - evt.s(0,2);
        double ds_new = std::pow(std::pow(dx, 2) + std::pow(dy, 2), 0.5);
//        if (ds_new < dS) {
//            dS = ds_new;
           // INFO("Got "<<dS<<" for "<<event.s(0,1)<<", "<<event.s(0,2)<< " at "<<i<<" out of "<<nEvents);
//        }
        dSs[i] = ds_new;
    }


    #pragma omp parallel for reduction(min:dS)
    for (int i=0;i<dSs.size();i++){
        double ds_new = dSs[i];
        if (ds_new < dS) dS = ds_new;
    }

    return dS;
}

int getBinFromRef(Event event, std::map<int, EventList> ref, int nBins, double fracNevents = 1){
    int binSign = 1;
    if (event.s(0,2) > event.s(0, 1)) binSign = -1;
    double dS = 100;
    int myBin =-1;
    for (int i=1;i<nBins + 1; i++){
        int bin = binSign * i;
        EventList evts = ref[bin];
        double ds_new = minS(event, evts, fracNevents);
        if (ds_new < dS){

            dS = ds_new;
            myBin = bin;
            //INFO("Got "<<ds_new <<" for bin "<<bin);
        }


    }

    return myBin;
}

std::vector<std::map<int, EventList> > binEventsFromRef(EventList sig, EventList tag, std::map<int, EventList> ref, int nBins, double fracNevents=1, bool DTag = false){
    std::vector<std::map<int, EventList> > v_binnedEvents;

    EventType sigType = sig.eventType();
    EventType tagType = tag.eventType();
    std::map<int, EventList> binnedSig;
    std::map<int, EventList> binnedTag;
    INFO("my EventType = "<<sigType<<" ref eventType = "<<ref[1].eventType());

    size_t nRefs=0;
    for (auto p:ref){ nRefs += p.second.size();}
    INFO("Have "<<nRefs<<" reference events");

    for (int i=1;i<nBins + 1; i++){
        int bin = i;
        binnedSig.insert(std::pair<int, EventList>({bin, EventList(sigType)}));
        binnedSig.insert(std::pair<int, EventList>({-bin, EventList(sigType)}));
        binnedTag.insert(std::pair<int, EventList>({bin, EventList(tagType)}));
        binnedTag.insert(std::pair<int, EventList>({-bin, EventList(tagType)}));
    }
    INFO("Binning "<<sig.size()<<" events");
    ProgressBar pb(60, trimmedString(__PRETTY_FUNCTION__));
    for (int i=0;i<sig.size();i++){
//        INFO(i<<" s01 = "<<sig[i].s(0,1));
        int myBin = getBinFromRef(sig[i], ref, nBins, fracNevents);
        binnedSig[myBin].push_back(sig[i]);

        if (DTag){
            int myBin_OT = getBinFromRef(tag[i], ref, nBins, fracNevents);
            binnedTag[myBin_OT].push_back(tag[i]);
            pb.print( double(i)/double(sig.size()), " of DT binning"  );

        }
        //INFO("I think the bin should be "<<myBin<< " for "<<sig[i].s(0,1) <<", "<<sig[i].s(0,2)<< "at "<<i<<" out of "<<sig.size());
        else{
            binnedTag[myBin].push_back(tag[i]);
            pb.print(double(i)/double(sig.size()), "of ST binning");
        }
    }
    pb.finish();

    v_binnedEvents.push_back(binnedSig);
    v_binnedEvents.push_back(binnedTag);
    return v_binnedEvents;

}

std::vector<std::map<std::pair<int, int>, EventList> > binDTEventsFromRef(EventList sig, EventList tag, std::map<int, EventList> ref, int nBins, double fracNevents=1, bool DTag = false){
    std::vector<std::map<std::pair<int, int>, EventList> > v_binnedEvents;
    std::map<std::pair<int, int>, EventList> binnedSig;
    std::map<std::pair<int, int>, EventList> binnedTag;
    std::vector<int> bins = {};
    for (int i=1;i<nBins+1;i++){
        bins.push_back(i);
        bins.push_back(-i);
    }
    EventType sigType = sig.eventType();
    EventType tagType = tag.eventType();

    ProgressBar pb(60, trimmedString(__PRETTY_FUNCTION__));
    for (int i=0;i<bins.size(); i++){
        for (int j=0;j<bins.size(); j++){
            int bin_sig = bins[i];
            int bin_tag = bins[j];
            std::pair<int, int> bin_pair({bin_sig, bin_tag});
            binnedSig.insert(std::pair<std::pair<int, int>, EventList>({bin_pair, EventList(sigType)}));
            binnedTag.insert(std::pair<std::pair<int, int>, EventList>({bin_pair, EventList(tagType)}));
        }
    }


    for (int i=0;i<sig.size();i++){
        int myBin = getBinFromRef(sig[i], ref, nBins, fracNevents);

       
        int myBin_OT = getBinFromRef(tag[i], ref, nBins, fracNevents);
        pb.print( double(i)/double(sig.size()), " of DT binning"  );

        std::pair<int, int> bin_pair({myBin, myBin_OT});
 
        binnedSig[bin_pair].push_back(sig[i]);
        binnedTag[bin_pair].push_back(tag[i]);

        }
 
    pb.finish();
    v_binnedEvents.push_back(binnedSig);
    v_binnedEvents.push_back(binnedTag);
    return v_binnedEvents;





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

std::vector<std::map<int, double> > cisi(CoherentSum A, CoherentSum C, std::map<int, EventList> map){
    std::vector<std::map<int, double> > r;
    std::map<int, double> c, s;
    r.push_back(c);
    r.push_back(s);
    for (auto p:map){
        int bin = p.first;
        EventList evts = p.second;
        complex_t z = cs_i (A, C, evts, evts.size());
        r[0].insert(std::pair<int, double>({bin, std::real(z)}));
        r[1].insert(std::pair<int, double>({bin, std::imag(z)}));
    }
    return r;

}



double fB(double f, double fbar, double x, double y, double c, double s){
    return f + fbar * (std::pow(x, 2) + std::pow(y, 2)) + 2 * std::pow(f*fbar, 0.5) *( c * x + s * y);
}

std::map<int, double> fBs(std::map<int, double> Fs, std::map<int, double> Fbars, double x, double y, std::map<int, double> c, std::map<int, double> s, int BSign){
    double sum = 0;
    std::map<int, double> m;
    for (auto p:Fs){
        int bin = p.first;
        double F = Fs[bin];
        double Fbar = Fbars[bin];
        double ci = c[bin];
        double si = s[bin];
        sum += fB(F, Fbar, x, y, ci, si);
    } 
    for (auto p:Fs){
        int bin = p.first;
        m.insert(std::pair<int, double>({bin, fB(Fs[bin], Fbars[bin], x, BSign* y, c[bin], s[bin])} ));
    }
    return m;
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


double getK(int i, CoherentSum A, std::map<int, EventList> m){
    double sum = 0;

    auto K = [A, &m](int i){
        double r = 0;
        for (auto evt : (*(&m))[i]){
            r += std::norm(A.getValNoCache(evt));
        }
        return r;
    };

    for (auto p : m){
        int bin = p.first;
        sum += K(bin);

    }
    return K(i)/sum;

}

complex_t getZ(int i, CoherentSum A, CoherentSum C, std::map<int, EventList> m){
    complex_t num = 0;
    double denA = 0;
    double denC = 0;
    EventList evts = m[i];
    for (auto evt: evts){
        num += A.getValNoCache(evt) * std::conj(C.getValNoCache(evt));
        denA += std::norm(A.getValNoCache(evt));
        denC += std::norm(C.getValNoCache(evt));
    }
    return num/std::pow(denA * denC, 0.5);
}

std::map<int, double> getKs(std::map<int, EventList> map, CoherentSum A){
    std::map<int, double> m;
    for (auto p:map){
        m.insert(std::pair<int, double> ({p.first, getK(p.first, A, map)}));
    }
    return m;
}


double getKGamma(int i, pCoherentSum A, std::map<int, EventList> m){
    double sum = 0;

    auto K = [A, &m](int i){
        double r = 0;
        for (auto evt : (*(&m))[i]){
            r += std::norm(A.getVal(evt));
        }
        return r;
    };

    for (auto p : m){
        int bin = p.first;
        sum += K(bin);

    }
    return K(i)/sum;

}

std::map<int, double> getKGammas(std::map<int, EventList> map, pCoherentSum A){
    std::map<int, double> m;
    for (auto p:map){
        m.insert(std::pair<int, double> ({p.first, getKGamma(p.first, A, map)}));
    }
    return m;
}


void testMI(MinuitParameterSet MPS, EventType eventType);
void GGSZ(MinuitParameterSet MPS, EventType eventType);
void testTim(MinuitParameterSet MPS, EventType eventType);
void genMI(MinuitParameterSet MPS, EventType eventType);
void fitMI(MinuitParameterSet MPS, EventType eventType);
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

  EventType eventType(pNames);
  if (NamedParameter<bool>("testTim", false)) testTim(MPS, eventType);
  if (NamedParameter<bool>("GGSZ", false))  GGSZ(MPS, eventType);
  if (NamedParameter<bool>("testMI", false))  testMI(MPS, eventType);
  if (NamedParameter<bool>("genMI", false))  genMI(MPS, eventType);
  if (NamedParameter<bool>("fitMI", false))  fitMI(MPS, eventType);
//
//
}
void testMI(MinuitParameterSet MPS, EventType eventType){
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

//        std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
        std::vector<int> bins;// = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
        for (int i=-nBins;i<0;i++) bins.push_back(i);
        for (int i=1;i<nBins+ 1;i++) bins.push_back(i);

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
        std::vector<int> bins;// = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
        for (int i=-nBins;i<0;i++) bins.push_back(i);
        for (int i=1;i<nBins+ 1;i++) bins.push_back(i);
//        std::vector<int> bins = {-8,-4,-2,-1, 1,2,4,8};
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
                _chi2 += std::pow((myBin - actualBin), 2);
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
        //Minimiser mini(fitMagAndBin, &MPS);
        Minimiser mini(chi2_binning, &MPS);
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
    
    if (NamedParameter<bool>("scanParam", true)){
        for (int idx=0;idx<bins.size();idx++){
            int bin = bins[idx];
            EventList evts = my_bins[idx];
            std::string paramName = NamedParameter<std::string>("scanParamName", "D0{K*(892)bar-{K0S0,pi-},pi+}_Re");
            int nSigma = NamedParameter<int>("nSigma", 2);
            int nSteps = NamedParameter<int>("nStepsMinForBin", 10);
            double p0 = MPS[paramName]->mean();
            double dp = MPS[paramName]->err();
            double pMin = p0 - nSigma * dp;
            double pMax = p0 + nSigma * dp;

            for (int i=0;i<nSteps;i++){
                
                auto evt = evts[0];
            
                MPS[paramName]->setCurrentFitVal(pMin + (pMax - pMin)*i/nSteps);
                A.prepare();
                C.prepare();
                int myBin = getBin(evt, nBins,  A, C);
                double dd = strongPhaseDiff(A, C, evt);
                INFO("Bin = "<<bin<<" myBin for "<<paramName<<" = "<<MPS[paramName]->mean()<<" is "<<myBin<<" dd = "<<dd);
            
            }
        }
    }

    if (NamedParameter<bool>("minForBin", false)){
    
        for (auto p : MPS){
            if (p->isFree()){
                double p0 = p->mean();
                double dp = p->err();

                int nSigma = NamedParameter<int>("nSigma", 2);
                int nSteps = NamedParameter<int>("nStepsMinForBin", 10);
                int nEvents = NamedParameter<int>("nEventsMinForBin", 100);
                double pMin = p0 - nSigma * dp;
                double pMax = p0 + nSigma * dp;
                for (int idx=0;idx<bins.size();idx++){
                    double bindiff=100;
                    int myNewBin = -1;
                    double pOpt = p0;
                    for (int i=0;i<nSteps;i++){
                        if (my_bins[idx].size()<nEvents) nEvents = my_bins[idx].size();
                        for (int j=0;j<nEvents;j++){
                        p->setCurrentFitVal(pMin + (pMax - pMin)*i/nSteps);
                        A.prepare();
                        C.prepare();
                        Event evt = my_bins[idx][j];
                        int myBin = getBin(evt,nBins, A, C);
                        double newBinDiff = std::abs(myBin - bins[idx]);
                        if (newBinDiff<bindiff) {
                            bindiff = newBinDiff;
                            pOpt = p->mean();
                            myNewBin = myBin;
                        }
                        
                        }
                    }
                    p->setCurrentFitVal(pOpt);
                    INFO("For bin "<<bins[idx]<<" I think "<<p->name()<<" should be "<<pOpt<<" whichi gives bin = "<<myNewBin<<" old value = "<<p0);

    

                }
               

            }
        }

        for (auto p : MPS){
            if (p->isFree()){
                cout<<p->name()<<" Free "<<p->mean()<<" "<<p->err()<<"\n";
            }
        }
    }

   if (NamedParameter<bool>("doPlot", true)){
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
    }



   
  
}





void GGSZ(MinuitParameterSet MPS, EventType eventType){
    int nBins = NamedParameter<int>("nBins", 8, "number of bins for MI");
  TRandom3 rndm;
  rndm.SetSeed( NamedParameter<size_t>("Seed",0) );
  gRandom = &rndm;
  int NIntBins = NamedParameter<int>("NIntBins", 1000);
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


    std::map<int, EventList> m_mc = binEvents(mc, nBins, A, C);


    std::vector<std::string>  KKstr = {"D0", "K+", "K-"};
    std::vector<std::string>  Kspi0str = {"D0", "K0S0", "pi0"};
    std::vector<std::string>  Kppimstr = {"D0", "K+", "pi-"};
    EventType KK(KKstr);
    EventType Kspi0(Kspi0str);
    EventType Kppim(Kppimstr);



    EventList sig_KK(eventType);
    EventList sig_Kppim(eventType);
    EventList sig_Kspipi(eventType);
    EventList tag_KK(KK);
    EventList tag_Kppim(Kppim);
    EventList tag_Kspipi(eventType);


	
    PhaseSpace phspSig(eventType,&rndm);
    PhaseSpace phspTag_KK(KK,&rndm);
    PhaseSpace phspTag_Kppim(Kppim,&rndm);
    pCorrelatedSum cs_KK(eventType, KK, MPS);
    pCorrelatedSum cs_Kppim(eventType, Kppim, MPS);
    pCorrelatedSum cs_Kspipi(eventType, eventType, MPS);
    INFO("Generating Events now!");
    size_t blockSize = NamedParameter<size_t>("blockSize", 10000);
    size_t nCorrEvents      = NamedParameter<size_t>     ("nCorrEvents"  , 10000, "Total number of events to generate" );
    GenerateCorrEvents( sig_KK, tag_KK, cs_KK, phspSig, phspTag_KK , nCorrEvents, blockSize, &rndm );
    GenerateCorrEvents( sig_Kppim, tag_Kppim, cs_Kppim, phspSig, phspTag_Kppim , nCorrEvents, blockSize, &rndm );
    GenerateCorrEvents( sig_Kspipi, tag_Kspipi, cs_Kspipi, phspSig, phspSig , nCorrEvents, blockSize, &rndm );

    pCoherentSum psiBplus( eventType, MPS, 1, true, true);
    pCoherentSum psiBminus( eventType, MPS, -1, true, false);
    EventList BpEvents(eventType);
    EventList BmEvents(eventType);

   
    size_t nBEvents      = NamedParameter<size_t>     ("nBEvents"  , 10000, "Total number of events to generate" );
    GenerateEvents(BpEvents, psiBplus, phspSig, nBEvents, blockSize, &rndm);
    GenerateEvents(BmEvents, psiBminus, phspSig, nBEvents, blockSize, &rndm);

    psiBplus.setEvents(BpEvents);
    psiBplus.setMC(mc);
    psiBminus.setEvents(BmEvents);
    psiBminus.setMC(mc);

    psiBplus.prepare();
    psiBminus.prepare();
 
    

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

    std::map<int, double> K = getKs(m_mc, A);
    std::map<int, double> Kbar = getKs(m_mc, C);
    std::map<int, double> KplusG = getKGammas(m_mc, psiBplus);
    std::map<int, double> KminusG = getKGammas(m_mc, psiBminus);

//    std::map<int, complex_t> cisi = cs(A, C, mc, nBins);
    
    std::vector<std::map<int, double> > _cs = cisi(A, C, m_mc);
    std::map<int, double> c = _cs[0];
    std::map<int, double> s = _cs[1];
    for (auto p : K){
        INFO("bin "<<p.first<<" K = "<<p.second<<" c = "<<c[p.first]<< " s = "<<s[p.first]);
    }

    auto chi2Bplus = [&K, &Kbar, &c, &s, MPS,nBEvents, &KplusG](){
        double chi2 = 0;
        double x = MPS["pCoherentSum::x+"]->mean();
        double y = MPS["pCoherentSum::y+"]->mean();
        std::map<int, double> FExp = fBs((*(&Kbar)) , (*(&K)),  x,  y, (*(&c)), (*(&s)), 1);
        for (auto p : FExp){
           // chi2 += std::pow(FExp[p.first] - (*(&KplusG))[p.first] , 2 );
            chi2 +=- 2 * (std::log(FExp[p.first] * nBEvents) * nBEvents * (*(&KplusG))[p.first]  - nBEvents * FExp[p.first]);

          // chi2 += -2 * std::log( FExp[p.first] );
        }
        return  chi2;
    };

    auto chi2Bminus = [&K, &Kbar, &c, &s, MPS,nBEvents, &KminusG](){
        double chi2 = 0;
        double x = MPS["pCoherentSum::x-"]->mean();
        double y = MPS["pCoherentSum::y-"]->mean();
        std::map<int, double> FExp = fBs((*(&K)), (*(&Kbar)) ,  x,  y, (*(&c)), (*(&s)), -1);
        for (auto p : FExp){
            //chi2 += std::pow(FExp[p.first] - (*(&KminusG))[p.first] , 2 );
            chi2 +=-2 *  (std::log(FExp[p.first] * nBEvents) * nBEvents * (*(&KminusG))[p.first]  - nBEvents * FExp[p.first]);
           //chi2 += -2 * std::log( FExp[p.first] );
        }
        return chi2;
    };

    auto chi2B = [&chi2Bplus, &chi2Bminus](){
        double chi2 = (*(&chi2Bplus))() + (*(&chi2Bminus))();
        return chi2;
    };

    INFO("chi2BP = "<<chi2Bplus());
    INFO("chi2BP = "<<chi2Bminus());

    Minimiser mini(chi2B, &MPS);
    mini.doFit();
    FitResult * fr = new FitResult(mini);
    fr->writeToFile(NamedParameter<std::string>("GGSZFitLog", "GGSZFit.log"));


    auto LL_pC = [&psiBplus, BpEvents, &psiBminus, BmEvents](){
        double ll = 0;
        //(&psiBplus)->prepare();
        //(&psiBminus)->prepare();
        double norm_plus = (&psiBplus)->norm();
        double norm_minus = (&psiBminus)->norm();
        for (auto evt : BpEvents){
            ll += -2* std::log(std::norm((&psiBplus)->getVal(evt))/norm_plus );
        }
        for (auto evt : BmEvents){
            ll += -2* std::log(std::norm((&psiBminus)->getVal(evt))/norm_minus );
        }
 
        return ll;
    };
    Minimiser miniLLBplus(LL_pC, &MPS);
    miniLLBplus.doFit();





}

void testTim(MinuitParameterSet MPS, EventType eventType){
 size_t hwt = std::thread::hardware_concurrency();
  size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
  //double luminosity   = NamedParameter<double>("Luminosity"  , 818.3       , "Luminosity to generate. Defaults to CLEO-c integrated luminosity.");
  //size_t nEvents      = NamedParameter<size_t>("nEvents"     , 0           , "Can also generate a fixed number of events per tag, if unspecified use the CLEO-c integrated luminosity.");
  size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
  //bool   poissonYield = NamedParameter<bool  >("PoissonYield", true        , "Flag to include Poisson fluctuations in expected yields (only if nEvents is not specified)");
  //bool   noQCFit      = NamedParameter<bool  >("noQCFit"     , false       , "Treat Signal and Tag as uncorrelated and fit the data individually");
  
    size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
    size_t nBins = NamedParameter<size_t>("nBins", 8);
    TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;


    std::vector<std::string>  KKstr = {"D0", "K+", "K-"};
    std::vector<std::string>  Kspi0str = {"D0", "K0S0", "pi0"};
    std::vector<std::string>  Kppimstr = {"D0", "K+", "pi-"};
    EventType KK(KKstr);
    EventType Kspi0(Kspi0str);
    EventType Kppim(Kppimstr);



    EventList sig_KK(eventType);
    EventList sig_Kppim(eventType);
    EventList sig_Kspipi(eventType);
    EventList tag_KK(KK);
    EventList tag_Kppim(Kppim);
    EventList tag_Kspipi(eventType);


	
    PhaseSpace phspSig(eventType,&rndm);
    PhaseSpace phspTag_KK(KK,&rndm);
    PhaseSpace phspTag_Kppim(Kppim,&rndm);
 //   pCorrelatedSum cs_KK(eventType, KK, MPS);


 




    std::string refBinningFile = NamedParameter<std::string>("refBinningFile", "ref_equal.root");

    TFile * fRefEqual = TFile::Open(refBinningFile.c_str());
    fRefEqual->cd();
    //std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
    std::vector<int> bins = {1,2,3,4,5,6,7,8, -1, -2, -3, -4, -5, -6, -7, -8};

    std::map<int, EventList> my_bins;
    for (int bin_idx=0;bin_idx<bins.size();bin_idx++){
        int bin = bins[bin_idx];
        INFO("bin = "<<bin);
        EventList my_bin_bin;
//       TTree * my_t_bin;
        if (bin>0) my_bin_bin = EventList(("ref_equal.root:myBin"+std::to_string(bin)).c_str(), eventType );
        if (bin<0) my_bin_bin = EventList(("ref_equal.root:myBinm"+std::to_string(std::abs(bin))).c_str(), eventType );
        my_bins.insert(std::pair<int, EventList>({bin, my_bin_bin}));

    }

    fRefEqual->Close();


    EventList mc =  Generator<>(eventType, &rndm).generate(NamedParameter<size_t>("nMC", 10000));
    CoherentSum A(eventType, MPS);
    CoherentSum C(eventType.conj(true), MPS);
    A.setEvents(mc);
    C.setEvents(mc);
    A.setMC(mc);
    C.setMC(mc);
    A.prepare();
    C.prepare();

    std::map<int, double> K;
    std::map<int, double> Kbar;
    std::map<int, complex_t> Z;
    for (int i=0;i<bins.size();i++){
            int bin = bins[i];
    //double getK(int i, CoherentSum A, std::map<int, EventList> m){
            double _K = getK(bin, A, my_bins);
            double _Kbar = getK(bin, C,  my_bins);
            complex_t _Z = getZ(bin, A, C, my_bins);
            K.insert(std::pair<int, double>({bin, _K}));
            Kbar.insert(std::pair<int, double>({bin, _Kbar}));
            Z.insert(std::pair<int, complex_t>({bin, _Z}));
    }


    if (NamedParameter<bool>("printKZ", false)){


        for (int i=0;i<bins.size();i++){
            int bin = bins[i];
    //double getK(int i, CoherentSum A, std::map<int, EventList> m){
//            double K = getK(bin, A, my_bins);
 //           double Kbar = getK(bin, C,  my_bins);

    //double getZ(int i, CoherentSum A, CoherentSum C, std::map<int, EventList> m){
//            complex_t Z = -getZ(bin, A, C, my_bins);
            INFO(bin<<" "<<std::real(Z[bin])<<" "<<std::imag(Z[bin]) <<" "<< K[bin] <<" "<<Kbar[bin] ) ;

        }


    }


    bool doCP = NamedParameter<bool>("doCP", true);
    if (doCP){
    pCorrelatedSum cs_KK(eventType, KK, MPS);
        GenerateCorrEvents( sig_KK, tag_KK, cs_KK, phspSig, phspTag_KK , NamedParameter<size_t>("nCorrEvents", 100000), blockSize, &rndm );

        EventList mcKK =  Generator<>(KK, &rndm).generate(NamedParameter<size_t>("nMC", 10000));
        std::vector<std::map<int, EventList> > binnedEvents = binEventsFromRef(sig_KK, tag_KK, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 0.1));
        std::map<int, double> NSig;
        std::map<int, double> NTag;
        std::map<int, double> YCP;
        std::map<int, double> pull;
        double sumNSig = 0;
        double sumNTag = 0;
        double sumYCP = 0;
        int CP = NamedParameter<int>("CP", 1);
        for (int i=0;i<bins.size();i++){
            int bin = bins[i];
            double NSig_i = binnedEvents[0][bin].size();
            double NTag_i = binnedEvents[1][bin].size();
            sumNSig += NSig_i;
            sumNTag += NTag_i;
            NSig.insert(std::pair<int, double>({bin, NSig_i} ));
            NTag.insert(std::pair<int, double>({bin, NTag_i} ));
//            double K = getK(bin, A, my_bins);
//            double Kbar = getK(bin, C,  my_bins);
//            complex_t Z = getZ(bin, A, C, my_bins);
            double _Y = K[bin] + Kbar[bin] - 2 * CP * std::pow(K[bin] * Kbar[bin],0.5) * std::real(Z[bin]);
            sumYCP += _Y;
            YCP.insert(std::pair<int,  double>({bin, _Y}));

 

        }
        for (int i=0;i<bins.size();i++){
            int bin = bins[i];
            double Y = YCP[bin]/sumYCP;
            double nSig = NSig[bin];
            double sSig = std::pow(nSig, 0.5);

            double p = (Y * sumNSig - nSig)/sSig;
            INFO(bin<<" "<<Y * sumNSig<<" "<<nSig<<" "<<p);
            pull.insert(std::pair<int, double>({bin, p}));


        }


        

        std::ofstream fPulls;
        fPulls.open(NamedParameter<std::string>("pullFileCP", "pullsCP.txt") ,std::ofstream::out | std::ofstream::app );
        for (int i=1;i<nBins+1;i++){
            fPulls<<pull[i]<<" "<<pull[-i]<<" ";
        }
        fPulls<<"\n";
        fPulls.close();
    }




    bool doKspipi = NamedParameter<bool>("doKspipi", true);
    if (doKspipi){
        pCorrelatedSum cs_Kspipi(eventType, eventType, MPS);
        EventList sig = sig_Kspipi;
        EventList tag = tag_Kspipi;
        PhaseSpace phspTag = phspSig;
        GenerateCorrEvents( sig, tag, cs_Kspipi, phspSig, phspTag , NamedParameter<size_t>("nCorrEvents", 100000), blockSize, &rndm );

        CoherentSum A(eventType, MPS);
        CoherentSum B(eventType.conj(true), MPS);
        CoherentSum C(eventType.conj(true), MPS);
        CoherentSum D(eventType, MPS);
        A.setEvents(sig);
        B.setEvents(tag);
        C.setEvents(sig);
        D.setEvents(tag);
        A.setMC(mc);
        B.setMC(mc);
        C.setMC(mc);
        D.setMC(mc);
        A.prepare();
        B.prepare();
        C.prepare();
        D.prepare();


        INFO("Binning DT");
        std::vector<std::map<std::pair<int, int>, EventList> > binnedEventsDT = binDTEventsFromRef(sig, tag, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 0.1));
        INFO("Binned DT");
        std::map<std::pair<int, int>, int> M;
        

        
        std::map<std::pair<int, int>, double> pullsKspipi;
        std::map<std::pair<int, int>, double> YKspipi;
        double sumYKspipi = 0;
        double sumM = 0;
        for (int i=0;i<bins.size();i++){
            for (int j=0;j<bins.size();j++){

                int bin_sig = bins[i];
                int bin_tag = bins[j];
                double K_sig = K[bin_sig];//getK(bin_sig, A, my_bins);
                double Kbar_sig = Kbar[bin_sig];//getK(bin_sig, C,  my_bins);
                double K_tag = K[bin_tag];//getK(bin_tag, A, my_bins);
                double Kbar_tag = Kbar[bin_tag];//getK(bin_tag, C,  my_bins);



        //double getZ(int i, CoherentSum A, CoherentSum C, std::map<int, EventList> m)
                complex_t Z_sig = Z[bin_sig];//getZ(bin_sig, A, C, my_bins);
                complex_t Z_tag = Z[bin_tag];//getZ(bin_tag, A, C, my_bins);

                //double Y = K * std::norm(b_Kspipi)  + Kbar * std::norm(d_Kspipi) - 2 * std::pow(K*Kbar, 0.5) * std::real(Z * z_Kspipi) ;
                double Y = K_sig * Kbar_tag + Kbar_sig * K_tag - 2 * std::pow(K_sig * K_tag * Kbar_sig * Kbar_tag, 0.5) * std::real(Z_sig * std::conj(Z_tag));
                //YKspipi.insert(std::pair<int, double>({bin , Y}));
                sumYKspipi += Y;
                std::pair<int, int> binPair({bin_sig, bin_tag});
                int M_ij = binnedEventsDT[0][binPair].size();
                M.insert(std::pair<std::pair<int, int> , double>({binPair, M_ij}));
                sumM += M_ij;
                YKspipi.insert(std::pair<std::pair<int, int> , double>({binPair, Y}));
                INFO(bin_sig<<", "<<bin_tag<<" "<<M_ij<<" "<<Y);
            }
        }
        INFO("Sum Mij = "<<sumM);
        for (auto p : YKspipi){

            int bin_sig = p.first.first;
            int bin_tag = p.first.second;

            double Y = p.second/sumYKspipi;
            double M_ij = M[p.first];
            double diff = Y * sumM - M_ij;
            double dM_ij = std::pow(M_ij, 0.5);

            double pull = -100;
            if (dM_ij != 0) pull = diff/dM_ij;
            INFO(bin_sig<<", "<<bin_tag<<" "<<Y * sumM<<" "<<M_ij<<" "<<diff<<" "<<pull);
            pullsKspipi.insert(std::pair<std::pair<int, int> , double>({p.first, pull}));

        }
        std::ofstream fPulls;
        fPulls.open(NamedParameter<std::string>("pullFileKspipi", "pullsKspipi.txt") ,std::ofstream::out | std::ofstream::app );
        for (auto p : pullsKspipi){
            fPulls<<p.second<<" ";
        }
        fPulls<<"\n";
        fPulls.close();
        }
    /*
    */
    bool doBPulls = NamedParameter<bool>("doBPulls", false);
    if (doBPulls){
        pCoherentSum psiBplus( eventType, MPS, 1, true, true);
        pCoherentSum psiBminus( eventType, MPS, -1, true, false);
        EventList BpEvents(eventType);
        EventList BmEvents(eventType);

    
        size_t nBEvents      = NamedParameter<size_t>     ("nBEvents"  , 1000, "Total number of events to generate" );
        GenerateEvents(BpEvents, psiBplus, phspSig, nBEvents, blockSize, &rndm);
        GenerateEvents(BmEvents, psiBminus, phspSig, nBEvents, blockSize, &rndm);

        std::vector<std::map<int, EventList> > binnedEventsBp = binEventsFromRef(BpEvents, BpEvents, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 0.1));
        std::vector<std::map<int, EventList> > binnedEventsBm = binEventsFromRef(BmEvents, BmEvents, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 0.1));

        std::map<int, double> Nplus;
        std::map<int, double> Nminus;
        double sumNplus = 0;
        double sumNminus = 0;
        for (int i=0;i<bins.size();i++){
            int bin = bins[i];
            double Nplus_i = binnedEventsBp[0][bin].size();
            double Nminus_i = binnedEventsBm[0][bin].size();
            sumNplus += Nplus_i;
            sumNminus += Nminus_i;
            Nplus.insert(std::pair<int, double>({bin, Nplus_i} ));
            Nminus.insert(std::pair<int, double>({bin, Nminus_i} ));

        }



        psiBplus.setEvents(BpEvents);
        psiBplus.setMC(mc);
        psiBminus.setEvents(BmEvents);
        psiBminus.setMC(mc);

        psiBplus.prepare();
        psiBminus.prepare();
        
        std::map<int, double> YBplus;
        std::map<int, double> YBminus;
        double sumYBplus = 0;
        double sumYBminus = 0;
        std::map<int, double> pullBplus;
        std::map<int, double> pullBminus;
        for (int i=0;i<bins.size();i++){

            int bin = bins[i];
            double K = getK(bin, A, my_bins);
            double Kbar = getK(bin, C,  my_bins);

    //double getZ(int i, CoherentSum A, CoherentSum C, std::map<int, EventList> m)
            complex_t Z = getZ(bin, A, C, my_bins);

            complex_t zBplus( MPS["pCoherentSum::x+"]->mean(), MPS["pCoherentSum::y+"]->mean() );
            complex_t zBminus( MPS["pCoherentSum::x-"]->mean(), MPS["pCoherentSum::y-"]->mean() );



            double _Yminus = K  + Kbar * std::norm(zBminus) + 2 * std::pow(K*Kbar, 0.5) * std::real(Z * std::conj(zBminus)) ;
            double _Yplus = Kbar  + K * std::norm(zBplus) + 2 * std::pow(K*Kbar, 0.5) * std::real(std::conj(Z) *std::conj(zBplus)) ;
            sumYBplus += _Yplus;
            sumYBminus += _Yminus;
            YBplus.insert(std::pair<int, double>({bin, _Yplus}));
            YBminus.insert(std::pair<int, double>({bin, _Yminus}));
        }
        for (int i=0;i<bins.size();i++){

            int bin = bins[i];
            double Yplus = YBplus[bin]/sumYBplus;
            double Yminus = YBminus[bin]/sumYBminus;
            double nplus = Nplus[bin];
            double nminus = Nminus[bin];
            double splus = std::pow(nplus, 0.5);
            double sminus = std::pow(nminus, 0.5);

            double pplus = (Yplus * sumNplus - nplus)/splus;
            double pminus = (Yminus * sumNminus - nminus)/sminus;
            INFO(bin<<" "<<pplus<<" "<<pminus<<" "<<Yplus<<" "<<nplus<<" "<<Yminus<<" "<<nminus);
            pullBplus.insert(std::pair<int, double>({bin, pplus}));
            pullBminus.insert(std::pair<int, double>({bin, pminus}));


        }
        std::ofstream fPullsBp;
        fPullsBp.open(NamedParameter<std::string>("pullFileBplus", "pullsBplus.txt") ,std::ofstream::out | std::ofstream::app );
        for (int i=1;i<nBins+1;i++){
            fPullsBp<<pullBplus[i]<<" "<<pullBplus[-i]<<" ";
        }
        fPullsBp<<"\n";
        fPullsBp.close();
        std::ofstream fPullsBm;
        fPullsBm.open(NamedParameter<std::string>("pullFileBminus", "pullsBminus.txt") ,std::ofstream::out | std::ofstream::app );
        for (int i=1;i<nBins+1;i++){
            fPullsBm<<pullBminus[i]<<" "<<pullBminus[-i]<<" ";
        }
        fPullsBm<<"\n";
        fPullsBm.close();

        }




//std::vector<std::map<int, EventList> > binEventsFromRef(EventList sig, EventList tag, std::map<int, EventList> ref, int nBins){

}



void genMI(MinuitParameterSet MPS, EventType eventType){

    size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
    
    TRandom3 rndm;
    rndm.SetSeed( NamedParameter<int>("Seed", 0) );
    gRandom = &rndm;
 

    size_t nBins(NamedParameter<size_t>("nBins", 8));
    std::vector<int> bins = {};
    for (int i=1;i<nBins+1;i++){
        bins.push_back(i);
        bins.push_back(-i);
    }


    size_t nD0(NamedParameter<size_t>("nD0", 5000));
    size_t nDbar0(NamedParameter<size_t>("nDbar0", 5000));
    size_t nCP(NamedParameter<size_t>("nCP", 5000));
    size_t nKspipi(NamedParameter<size_t>("nKspipi", 1000));
    size_t nB(NamedParameter<size_t>("nB", 2000));
    
 
    std::vector<std::string>  KKstr = {"D0", "K+", "K-"};
    std::vector<std::string>  Kspi0str = {"D0", "K0S0", "pi0"};
    std::vector<std::string>  Kppimstr = {"D0", "K+", "pi-"};
    std::vector<std::string>  Kmpipstr = {"D0", "K-", "pi+"};
    EventType KK(KKstr);
    EventType Kspi0(Kspi0str);
    EventType Kppim(Kppimstr);
    EventType Kmpip(Kmpipstr);

    EventList sig_KK(eventType);
    EventList sig_Kspi0(eventType);
    EventList sig_Kppim(eventType);
    EventList sig_Kmpip(eventType);
    EventList sig_Kspipi(eventType);
    EventList tag_KK(KK);
    EventList tag_Kspi0(Kspi0);
    EventList tag_Kppim(Kppim);
    EventList tag_Kmpip(Kmpip);
    EventList tag_Kspipi(eventType);

    EventList sig_Bp(eventType);
    EventList sig_Bm(eventType);
 


    PhaseSpace phspSig(eventType,&rndm);
    PhaseSpace phspTag_KK(KK,&rndm);
    PhaseSpace phspTag_Kspi0(Kspi0,&rndm);
    PhaseSpace phspTag_Kppim(Kppim,&rndm);
    PhaseSpace phspTag_Kmpip(Kmpip,&rndm);


    INFO("Generating correlated events");

    pCorrelatedSum cs_KK(eventType, KK, MPS);
    pCorrelatedSum cs_Kmpip(eventType, Kmpip, MPS);
    pCorrelatedSum cs_Kppim(eventType, Kppim, MPS);
    pCorrelatedSum cs_Kspi0(eventType, Kspi0, MPS);
    pCorrelatedSum cs_Kspipi(eventType, eventType, MPS);

    std::string refBinningFile = NamedParameter<std::string>("refBinningFile", "ref_equal.root");

    TFile * fRefEqual = TFile::Open(refBinningFile.c_str());
    fRefEqual->cd();
    //std::vector<int> bins = {-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8};
    //std::vector<int> bins = {1,2,3,4,5,6,7,8, -1, -2, -3, -4, -5, -6, -7, -8};

    std::map<int, EventList> my_bins;
    for (int bin_idx=0;bin_idx<bins.size();bin_idx++){
        int bin = bins[bin_idx];
        INFO("bin = "<<bin);
        EventList my_bin_bin;
//       TTree * my_t_bin;
        if (bin>0) my_bin_bin = EventList(("ref_equal.root:myBin"+std::to_string(bin)).c_str(), eventType );
        if (bin<0) my_bin_bin = EventList(("ref_equal.root:myBinm"+std::to_string(std::abs(bin))).c_str(), eventType );
        my_bins.insert(std::pair<int, EventList>({bin, my_bin_bin}));

    }

    fRefEqual->Close();



    std::string output = NamedParameter<std::string>("outputMI", "binnedEventsForMI.root");
    TFile * f = TFile::Open(output.c_str(), "UPDATE");
    f->cd();
    std::string treeName_KK = "KK";
    std::string treeName_Kspi0 = "Kspi0";
    std::string treeName_Kppim = "Kppim";
    std::string treeName_Kmpip = "Kmpip";
    std::string treeName_Bp = "Bp";
    std::string treeName_Bm = "Bm";
 

    
    EventList mc =  Generator<>(eventType, &rndm).generate(NamedParameter<size_t>("nMC", 10000));
    CoherentSum A(eventType, MPS);
    CoherentSum C(eventType.conj(true), MPS);
    A.setEvents(mc);
    C.setEvents(mc);
    A.setMC(mc);
    C.setMC(mc);
    A.prepare();
    C.prepare();


    INFO("Getting K_i, Z_i from the model");
    std::map<int, double> K;
    std::map<int, double> Kbar;
    std::map<int, complex_t> Z;
    for (int i=0;i<bins.size();i++){
            int bin = bins[i];
            double _K = getK(bin, A, my_bins);
            double _Kbar = getK(bin, C,  my_bins);
            complex_t _Z = getZ(bin, A, C, my_bins);
            K.insert(std::pair<int, double>({bin, _K}));
            Kbar.insert(std::pair<int, double>({bin, _Kbar}));
            Z.insert(std::pair<int, complex_t>({bin, _Z}));
    }
    
    INFO("Generating KK events");
    GenerateCorrEvents( sig_KK, tag_KK, cs_KK, phspSig, phspTag_KK , nCP, blockSize, &rndm );
    INFO("Binning KK << nEvents = "<<sig_KK.size());
    std::vector<std::map<int, EventList> > binned_KK = binEventsFromRef(sig_KK, tag_KK, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));

    for (int i=0;i<bins.size();i++){
        int bin = bins[i];
        INFO(bin<<" "<<binned_KK[0][bin].size());
   

    }


 
    INFO("Generating Kspipi events");
    GenerateCorrEvents( sig_Kspipi, tag_Kspipi, cs_Kspipi, phspSig, phspSig , nKspipi, blockSize, &rndm );
    INFO("Binning Kspipi");
    std::vector<std::map<std::pair<int, int>, EventList> > binned_Kspipi = binDTEventsFromRef(sig_Kspipi, tag_Kspipi, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));


    std::string treeName_Kspipi = "Kspipi";
    for (auto p: binned_Kspipi[0]){
        std::pair<int, int> bin_pair = p.first;
        int bin_sig = bin_pair.first;
        int bin_tag = bin_pair.second;
        EventList events_sig = binned_Kspipi[0][bin_pair];
        EventList events_tag = binned_Kspipi[1][bin_pair];
        INFO("n("<<bin_sig<<", "<<bin_tag<<") = "<<events_sig.size()<<" "<<events_sig.eventType()<<" "<<events_tag.eventType());
        
        std::string treeName_Kspipi_sig_ij = treeName_Kspipi + "_sig_" + std::to_string(bin_sig) + "_" + std::to_string(bin_tag);
        std::string treeName_Kspipi_tag_ij = treeName_Kspipi + "_tag_" + std::to_string(bin_sig) + "_" + std::to_string(bin_tag);
        if (events_sig.size()!= 0){
            events_sig.tree(treeName_Kspipi_sig_ij.c_str())->Write(treeName_Kspipi_sig_ij.c_str());        
            events_tag.tree(treeName_Kspipi_tag_ij.c_str())->Write(treeName_Kspipi_tag_ij.c_str());        
        }
    }




    INFO("Generating Kspi0 events");
    GenerateCorrEvents( sig_Kspi0, tag_Kspi0, cs_Kspi0, phspSig, phspTag_Kspi0 , nCP, blockSize, &rndm );
    INFO("Generating Kppim events");
    GenerateCorrEvents( sig_Kppim, tag_Kppim, cs_Kppim, phspSig, phspTag_Kppim , nD0, blockSize, &rndm );
    INFO("Generating Kmpip events");
    GenerateCorrEvents( sig_Kmpip, tag_Kmpip, cs_Kmpip, phspSig, phspTag_Kmpip , nDbar0, blockSize, &rndm );
    
    
    INFO("Generating B+- events");
    pCoherentSum psiBplus( eventType, MPS, 1, true, true);
    pCoherentSum psiBminus( eventType, MPS, -1, true, false);

    GenerateEvents(sig_Bp, psiBplus, phspSig, nB, blockSize, &rndm);
    GenerateEvents(sig_Bm, psiBminus, phspSig, nB, blockSize, &rndm);


    INFO("Binning Kspi0");
    std::vector<std::map<int, EventList> > binned_Kspi0 = binEventsFromRef(sig_Kspi0, tag_Kspi0, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));
    INFO("Binning Kppim");
    std::vector<std::map<int, EventList> > binned_Kppim = binEventsFromRef(sig_Kppim, tag_Kppim, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));
    INFO("Binning Kmpip");
    std::vector<std::map<int, EventList> > binned_Kmpip = binEventsFromRef(sig_Kmpip, tag_Kmpip, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));

    INFO("Binning B+");
    std::vector<std::map<int, EventList> > binned_Bp = binEventsFromRef(sig_Bp, sig_Bp, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));
    INFO("Binning B-");
    std::vector<std::map<int, EventList> > binned_Bm = binEventsFromRef(sig_Bm, sig_Bm, my_bins, nBins, NamedParameter<double>("nEventsRefFrac", 1));

   
    for (int i=0;i<bins.size();i++){
        int bin = bins[i];
        std::string treeName_KK_sig_i = treeName_KK + "_sig_" +  std::to_string(bin);
        std::string treeName_Kppim_sig_i = treeName_Kppim + "_sig_" +  std::to_string(bin);
        std::string treeName_Kmpip_sig_i = treeName_Kmpip + "_sig_" +  std::to_string(bin);
        std::string treeName_Kspi0_sig_i = treeName_Kspi0 + "_sig_" +  std::to_string(bin);

        std::string treeName_KK_tag_i = treeName_KK + "_tag_" +  std::to_string(bin);
        std::string treeName_Kppim_tag_i = treeName_Kppim + "_tag_" +  std::to_string(bin);
        std::string treeName_Kmpip_tag_i = treeName_Kmpip + "_tag_" +  std::to_string(bin);
        std::string treeName_Kspi0_tag_i = treeName_Kspi0 + "_tag_" +  std::to_string(bin);


        std::string treeName_Bp_i = treeName_Bp + "_" +  std::to_string(bin);
        std::string treeName_Bm_i = treeName_Bm + "_" +  std::to_string(bin);

        INFO(bin<<" "<<binned_KK[0][bin].size()<<" "<<binned_Kspi0[0][bin].size());


        if (binned_KK[0][bin].size() != 0) binned_KK[0][bin].tree(treeName_KK_sig_i.c_str())->Write(treeName_KK_sig_i.c_str());
        if (binned_Kppim[0][bin].size() != 0) binned_Kppim[0][bin].tree(treeName_Kppim_sig_i.c_str())->Write(treeName_Kppim_sig_i.c_str());
        if (binned_Kmpip[0][bin].size() != 0) binned_Kmpip[0][bin].tree(treeName_Kmpip_sig_i.c_str())->Write(treeName_Kmpip_sig_i.c_str());
        if (binned_Kspi0[0][bin].size() != 0) binned_Kspi0[0][bin].tree(treeName_Kspi0_sig_i.c_str())->Write(treeName_Kspi0_sig_i.c_str());
        if (binned_KK[1][bin].size() != 0) binned_KK[1][bin].tree(treeName_KK_tag_i.c_str())->Write(treeName_KK_tag_i.c_str());
        if (binned_Kppim[1][bin].size() != 0) binned_Kppim[1][bin].tree(treeName_Kppim_tag_i.c_str())->Write(treeName_Kppim_tag_i.c_str());
        if (binned_Kmpip[1][bin].size() != 0) binned_Kmpip[1][bin].tree(treeName_Kmpip_tag_i.c_str())->Write(treeName_Kmpip_tag_i.c_str());
        if (binned_Kspi0[1][bin].size() != 0) binned_Kspi0[1][bin].tree(treeName_Kspi0_tag_i.c_str())->Write(treeName_Kspi0_tag_i.c_str());
        if (binned_Bp[0][bin].size() != 0) binned_Bp[0][bin].tree(treeName_Bp_i.c_str())->Write(treeName_Bp_i.c_str());
        if (binned_Bm[0][bin].size() != 0) binned_Bm[0][bin].tree(treeName_Bm_i.c_str())->Write(treeName_Bm_i.c_str());
    }
    f->Close();


}


void fitMI(MinuitParameterSet MPS, EventType eventType){
    TRandom3 rndm;
    rndm.SetSeed( NamedParameter<int>("Seed", 0) );
    gRandom = &rndm;
 

    size_t nBins(NamedParameter<size_t>("nBins", 8));
    std::vector<int> bins = {};
    for (int i=1;i<nBins+1;i++){
        bins.push_back(i);
        bins.push_back(-i);
    }


    std::string input = NamedParameter<std::string>("DataSample", "binnedEventsForMI.root");
    TFile * f = TFile::Open(input.c_str());
    std::map<int, EventList> binned_KK_sig;
    std::map<int, EventList> binned_Kppim_sig;
    std::map<int, EventList> binned_Kmpip_sig;
    std::map<int, EventList> binned_Kspi0_sig;
    std::map<std::pair<int, int>, EventList> binned_Kspipi_sig;
    std::map<int, EventList> binned_KK_tag;
    std::map<int, EventList> binned_Kppim_tag;
    std::map<int, EventList> binned_Kmpip_tag;
    std::map<int, EventList> binned_Kspi0_tag;
    std::map<std::pair<int, int>, EventList> binned_Kspipi_tag;
    
    
    std::map<int, EventList> binned_Bp;
    std::map<int, EventList> binned_Bm;
    std::string treeName_KK = "KK";
    std::string treeName_Kspi0 = "Kspi0";
    std::string treeName_Kppim = "Kppim";
    std::string treeName_Kmpip = "Kmpip";
    std::string treeName_Kspipi = "Kspipi";
    std::string treeName_Bp = "Bp";
    std::string treeName_Bm = "Bm";
 
    std::vector<std::string>  KKstr = {"D0", "K+", "K-"};
    std::vector<std::string>  Kspi0str = {"D0", "K0S0", "pi0"};
    std::vector<std::string>  Kppimstr = {"D0", "K+", "pi-"};
    std::vector<std::string>  Kmpipstr = {"D0", "K-", "pi+"};
    EventType KK(KKstr);
    EventType Kspi0(Kspi0str);
    EventType Kppim(Kppimstr);
    EventType Kmpip(Kmpipstr);

    EventList sig_KK(eventType);
    EventList sig_Kspi0(eventType);
    EventList sig_Kppim(eventType);
    EventList sig_Kmpip(eventType);
    EventList sig_Kspipi(eventType);
    EventList tag_KK(KK);
    EventList tag_Kspi0(Kspi0);
    EventList tag_Kppim(Kppim);
    EventList tag_Kmpip(Kmpip);
    EventList tag_Kspipi(eventType);

    EventList sig_Bp(eventType);
    EventList sig_Bm(eventType);
    

    for (int i=0;i<bins.size();i++){
        int bin = bins[i];
        std::string treeName_KK_sig_i = treeName_KK + "_sig_" +  std::to_string(bin);
        std::string treeName_Kppim_sig_i = treeName_Kppim + "_sig_" +  std::to_string(bin);
        std::string treeName_Kmpip_sig_i = treeName_Kmpip + "_sig_" +  std::to_string(bin);
        std::string treeName_Kspi0_sig_i = treeName_Kspi0 + "_sig_" +  std::to_string(bin);

        std::string treeName_KK_tag_i = treeName_KK + "_tag_" +  std::to_string(bin);
        std::string treeName_Kppim_tag_i = treeName_Kppim + "_tag_" +  std::to_string(bin);
        std::string treeName_Kmpip_tag_i = treeName_Kmpip + "_tag_" +  std::to_string(bin);
        std::string treeName_Kspi0_tag_i = treeName_Kspi0 + "_tag_" +  std::to_string(bin);


        std::string treeName_Bp_i = treeName_Bp + "_" +  std::to_string(bin);
        std::string treeName_Bm_i = treeName_Bm + "_" +  std::to_string(bin);

        INFO("Getting binned KK events "<<bin);
        if (f->Get(treeName_KK_sig_i.c_str()) != nullptr){
            TTree * tree_KK_sig_i = (TTree*)f->Get(treeName_KK_sig_i.c_str());
            TTree * tree_KK_tag_i = (TTree*)f->Get(treeName_KK_tag_i.c_str());
            sig_KK = EventList(tree_KK_sig_i, eventType);
            tag_KK = EventList(tree_KK_tag_i, KK);
        }
        else{
            sig_KK=EventList(eventType);
            tag_KK=EventList(KK);
        }



        INFO("Getting binned Kppim events "<<bin);
        if (f->Get(treeName_Kppim_sig_i.c_str()) != nullptr){
            TTree * tree_Kppim_sig_i = (TTree*)f->Get(treeName_Kppim_sig_i.c_str());
            TTree * tree_Kppim_tag_i = (TTree*)f->Get(treeName_Kppim_tag_i.c_str());
            sig_Kppim = EventList(tree_Kppim_sig_i, eventType);
            tag_Kppim = EventList(tree_Kppim_tag_i, Kppim);
        }
        else{
            sig_Kppim=EventList(eventType);
            tag_Kppim=EventList(Kppim);
        }

        INFO("Getting binned Kmpip events "<<bin);
        if (f->Get(treeName_Kmpip_sig_i.c_str()) != nullptr){
            TTree * tree_Kmpip_sig_i = (TTree*)f->Get(treeName_Kmpip_sig_i.c_str());
            TTree * tree_Kmpip_tag_i = (TTree*)f->Get(treeName_Kmpip_tag_i.c_str());
            sig_Kmpip = EventList(tree_Kmpip_sig_i, eventType);
            tag_Kmpip = EventList(tree_Kmpip_tag_i, Kmpip);
        }
        else{
            sig_Kmpip=EventList(eventType);
            tag_Kmpip=EventList(Kmpip);
            }

        INFO("Getting binned Kspi0 events "<<bin);
        if (f->Get(treeName_Kspi0_sig_i.c_str()) != nullptr){
            TTree * tree_Kspi0_sig_i = (TTree*)f->Get(treeName_Kspi0_sig_i.c_str());
            TTree * tree_Kspi0_tag_i = (TTree*)f->Get(treeName_Kspi0_tag_i.c_str());
            sig_Kspi0 = EventList(tree_Kspi0_sig_i, eventType);
            tag_Kspi0 = EventList(tree_Kspi0_tag_i, Kspi0);
        }
        else{
            sig_Kspi0=EventList(eventType);
            tag_Kspi0=EventList(Kspi0);
        }

        INFO("Getting binned B+ events "<<bin);
        if (f->Get(treeName_Bp_i.c_str()) != nullptr){
            TTree * tree_Bp_i = (TTree*)f->Get(treeName_Bp_i.c_str());
            sig_Bp = EventList(tree_Bp_i, eventType);
        }
        else{
            sig_Bp=EventList(eventType);

        }

        INFO("Getting binned B- events "<<bin);
        if (f->Get(treeName_Bm_i.c_str()) != nullptr){
            TTree * tree_Bm_i = (TTree*)f->Get(treeName_Bm_i.c_str());
            sig_Bm = EventList(tree_Bm_i, eventType);
        }
        else{
            sig_Bm=EventList(eventType);

        }




        binned_KK_sig.insert(std::pair<int, EventList>({bin, sig_KK}));
        binned_Kppim_sig.insert(std::pair<int, EventList>({bin, sig_Kppim}));
        binned_Kmpip_sig.insert(std::pair<int, EventList>({bin, sig_Kmpip}));
        binned_Kspi0_sig.insert(std::pair<int, EventList>({bin, sig_Kspi0}));
        binned_KK_tag.insert(std::pair<int, EventList>({bin, tag_KK}));
        binned_Kppim_tag.insert(std::pair<int, EventList>({bin, tag_Kppim}));
        binned_Kmpip_tag.insert(std::pair<int, EventList>({bin, tag_Kmpip}));
        binned_Kspi0_tag.insert(std::pair<int, EventList>({bin, tag_Kspi0}));

        binned_Bp.insert(std::pair<int, EventList>({bin, sig_Bp}));
        binned_Bm.insert(std::pair<int, EventList>({bin, sig_Bm}));


        for (int j=0;j<bins.size();j++){
            int bin_tag = bins[j];
            INFO("Getting binned Kspipi events "<<bin<<", "<<bin_tag);
            std::pair<int, int> binPair({bin, bin_tag});
            std::string treeName_Kspipi_sig_ij = treeName_Kspipi + "_sig_" + std::to_string(bin) + "_" + std::to_string(bin_tag);
            std::string treeName_Kspipi_tag_ij = treeName_Kspipi + "_tag_" + std::to_string(bin) + "_" + std::to_string(bin_tag);

            if (f->Get(treeName_Kspipi_sig_ij.c_str()) != nullptr){
                TTree * tree_Kspipi_sig_ij = (TTree*)f->Get(treeName_Kspipi_sig_ij.c_str());
                TTree * tree_Kspipi_tag_ij = (TTree*)f->Get(treeName_Kspipi_tag_ij.c_str());

                sig_Kspipi = EventList(tree_Kspipi_sig_ij, eventType);
                tag_Kspipi = EventList(tree_Kspipi_tag_ij, eventType);
            }
            else{
                sig_Kspipi = EventList(eventType);
                tag_Kspipi = EventList(eventType);
            }
            binned_Kspipi_sig.insert(std::pair<std::pair<int, int>, EventList>({binPair, sig_Kspipi})  );
            binned_Kspipi_tag.insert(std::pair<std::pair<int, int>, EventList>({binPair, tag_Kspipi})  );

        }



    }
    

    f->Close();

    std::map<int, double> K;
    std::map<int, double> Kbar;
    double sumK=0;
    double sumKbar = 0;
    std::map<int, double> F_KK;
    double sum_KK = 0;
    std::map<int, double> F_Kspi0;
    double sum_Kspi0 = 0;

    std::map<int, double> F_Bp;
    std::map<int, double> F_Bm;
    double sum_Bp = 0;
    double sum_Bm = 0;


    for (auto p : binned_Kppim_sig){
        double nD0_i = binned_Kppim_sig[p.first].size();
        double nDbar0_i = binned_Kppim_sig[p.first].size();
        sumK += nD0_i;
        sumKbar += nDbar0_i;
        K.insert(std::pair<int, double>({p.first, nD0_i})  );
        Kbar.insert(std::pair<int, double>({p.first, nDbar0_i})  );
        double nKK_i = binned_KK_sig[p.first].size();
        sum_KK += nKK_i;
        double nKspi0_i = binned_Kspi0_sig[p.first].size();
        sum_Kspi0 += nKspi0_i;
        F_KK.insert(std::pair<int, double>({p.first, nKK_i}));
        F_Kspi0.insert(std::pair<int, double>({p.first, nKspi0_i}));

        double nBp_i = binned_Bp[p.first].size();
        double nBm_i = binned_Bm[p.first].size();
        F_Bp.insert(std::pair<int, double>({p.first, nBp_i}));
        F_Bm.insert(std::pair<int, double>({p.first, nBm_i}));
        sum_Bp += nBp_i;
        sum_Bm += nBm_i;

    }
    for (auto p : K){
        K[p.first] = K[p.first]/sumK;
        Kbar[p.first] = Kbar[p.first]/sumKbar;
        F_KK[p.first] = F_KK[p.first]/sum_KK;
        F_Kspi0[p.first] = F_Kspi0[p.first]/sum_Kspi0;
        F_Bp[p.first] = F_Bp[p.first]/sum_Bp;
        F_Bm[p.first] = F_Bm[p.first]/sum_Bm;
    }

    std::map<std::pair<int, int>, double> F_Kspipi;
    double sum_Kspipi = 0;
    for (auto p : binned_Kspipi_sig){
        double n_Kspipi_ij = p.second.size();
        sum_Kspipi += n_Kspipi_ij;
        F_Kspipi.insert(std::pair<std::pair<int, int>, double>(p.first, n_Kspipi_ij));
    }
    for (auto p : binned_Kspipi_sig){
        F_Kspipi[p.first] = F_Kspipi[p.first]/sum_Kspipi;
    }

    auto Y_CP = [MPS](int bin, int CP, std::map<int, double> K, std::map<int, double> Kbar){
        double sumYCP = 0;

        for (auto p : K ){
            sumYCP += K[p.first] + Kbar[p.first] - 2 * std::pow(K[p.first]  * Kbar[p.first], 0.5  ) * CP * MPS["c_" + std::to_string(std::abs(p.first))]->mean();
        }
        
        return (K[bin] + Kbar[bin] - 2 * CP * std::pow(K[bin]* Kbar[bin], 0.5) * MPS["c_" + std::to_string(std::abs(bin))]->mean())/sumYCP;

    };

    auto Y_Kspipi = [MPS](int bin_sig, int bin_tag, std::map<int, double> K, std::map<int, double> Kbar){
        double sumY = 0;
        for (auto p1 : K){
            for (auto p2 : K){
                sumY += K[p1.first] * Kbar[p2.first] + Kbar[p1.first] * K[p2.first] 
                - 2 * std::pow( K[p1.first] * Kbar[p2.first] * Kbar[p1.first] * K[p2.first]  , 0.5 ) * (
                    MPS["c_" + std::to_string(std::abs(p1.first))]->mean() * MPS["c_" + std::to_string(std::abs(p2.first))]->mean()
                    +  p1.first * p2.first/(std::abs(p1.first) * std::abs(p2.first)) * MPS["s_" + std::to_string(std::abs(p1.first))]->mean() * MPS["s_" + std::to_string(std::abs(p2.first))]->mean()
                );

            }
        }
        return (K[bin_sig] * Kbar[bin_tag] + Kbar[bin_sig] * K[bin_tag] 
            - 2 * std::pow( K[bin_sig] * Kbar[bin_tag] * Kbar[bin_sig] * K[bin_tag]  , 0.5 ) * (
                MPS["c_" + std::to_string(std::abs(bin_sig))]->mean() * MPS["c_" + std::to_string(std::abs(bin_tag))]->mean()
                +  bin_sig * bin_tag/(std::abs(bin_sig) * std::abs(bin_tag)) * MPS["s_" + std::to_string(std::abs(bin_sig))]->mean() * MPS["s_" + std::to_string(std::abs(bin_tag))]->mean()
            ))/sumY;

    };
    auto Y_Bp = [MPS](int bin, std::map<int, double> K, std::map<int, double> Kbar, std::map<int, complex_t> Z){
        double sumY = 0;



        complex_t zB(MPS["pCoherentSum::x+"]->mean() ,-MPS["pCoherentSum::y+"]->mean());


        for (auto p : K ){
            //complex_t zC(MPS["c_" + std::to_string(std::abs(p.first)) ]->mean(), p.first/std::abs(p.first)  * MPS["s_" + std::to_string(std::abs(p.first)) ]->mean())

            complex_t zC =Z[p.first] ;

            double _Y = Kbar[p.first] + K[p.first] *(std::norm(zB)) + 2 * std::pow(K[p.first]  * Kbar[p.first], 0.5  ) * std::real(zC * std::conj(zB));
            if (_Y < 0){
                INFO("Y  = "<<_Y);
                _Y =0 ;

            }

            else{
            sumY += _Y;
            }
        }
        

//        complex_t zC(MPS["c_" + std::to_string(std::abs(bin)) ]->mean(), bin/std::abs(bin)  * MPS["s_" + std::to_string(std::abs(bin)) ]->mean());
        complex_t zC = Z[bin];
        double Y =  (Kbar[bin] + K[bin] *(std::norm(zB)) + 2 * std::pow(K[bin]  * Kbar[bin], 0.5  ) * std::real(zC * std::conj(zB)))/sumY;
        if (Y<0){
            INFO("Y = "<<Y);
            Y = 0;
        }
        return Y;

    };

    auto Y_Bm = [MPS](int bin, std::map<int, double> K, std::map<int, double> Kbar, std::map<int, complex_t> Z){
        double sumY = 0;
        //A + zC ==> A2 + C2 z2 + 2  (sqrt(A2C2) real(c + is) z*)

        complex_t zB(MPS["pCoherentSum::x-"]->mean() ,MPS["pCoherentSum::y-"]->mean());
//       complex_t zB(x, y);


        for (auto p : K ){
            //complex_t zC(MPS["c_" + std::to_string(std::abs(p.first)) ]->mean(), p.first/std::abs(p.first)  * MPS["s_" + std::to_string(std::abs(p.first)) ]->mean());

            complex_t zC =Z[p.first] ;
            double _Y = K[p.first] + Kbar[p.first] *(std::norm(zB)) + 2 * std::pow(K[p.first]  * Kbar[p.first], 0.5  ) * std::real(zC * std::conj(zB));
            
            if (_Y < 0){
                INFO("Y  = "<<_Y);
                _Y =0 ;

            }

            else{
            sumY += _Y;
            }
        }
        

        //complex_t zC(MPS["c_" + std::to_string(std::abs(bin)) ]->mean(), bin/std::abs(bin)  * MPS["s_" + std::to_string(std::abs(bin)) ]->mean());

        complex_t zC = Z[bin];
        double Y =  (K[bin] + Kbar[bin] *(std::norm(zB)) + 2 * std::pow(K[bin]  * Kbar[bin], 0.5  ) * std::real(zC * std::conj(zB)))/sumY;
        if (Y<0){
            INFO("Y = "<<Y);
            Y = 0;
        }
        return Y;




    };





    std::string refBinningFile = NamedParameter<std::string>("refBinningFile", "ref_equal.root");

    TFile * fRefEqual = TFile::Open(refBinningFile.c_str());
    fRefEqual->cd();
   
    

    std::map<int, EventList> my_bins;
    for (int bin_idx=0;bin_idx<bins.size();bin_idx++){
        int bin = bins[bin_idx];
        INFO("bin = "<<bin);
        EventList my_bin_bin;
//       TTree * my_t_bin;
        if (bin>0) my_bin_bin = EventList(("ref_equal.root:myBin"+std::to_string(bin)).c_str(), eventType );
        if (bin<0) my_bin_bin = EventList(("ref_equal.root:myBinm"+std::to_string(std::abs(bin))).c_str(), eventType );
        my_bins.insert(std::pair<int, EventList>({bin, my_bin_bin}));

    }

    fRefEqual->Close();

        
    EventList mc =  Generator<>(eventType, &rndm).generate(NamedParameter<size_t>("nMC", 10000));
    CoherentSum A(eventType, MPS);
    CoherentSum C(eventType.conj(true), MPS);
    A.setEvents(mc);
    C.setEvents(mc);
    A.setMC(mc);
    C.setMC(mc);
    A.prepare();
    C.prepare();


    INFO("Getting K_i, Z_i from the model");
    std::map<int, double> K0;
    std::map<int, double> Kbar0;
    std::map<int, complex_t> Z0;
    for (int i=0;i<bins.size();i++){
            int bin = bins[i];
            double _K = getK(bin, A, my_bins);
            double _Kbar = getK(bin, C,  my_bins);
            complex_t _Z = getZ(bin, A, C, my_bins);
            K0.insert(std::pair<int, double>({bin, _K}));
            Kbar0.insert(std::pair<int, double>({bin, _Kbar}));
            Z0.insert(std::pair<int, complex_t>({bin, _Z}));

    }

    for (int i=1;i<nBins + 1;i++){
        MPS["c_"+std::to_string(i)]->setCurrentFitVal(std::real(Z0[i]));
        MPS["s_"+std::to_string(i)]->setCurrentFitVal(std::imag(Z0[i]));
    }


    for (auto p : F_KK){
        double N_KK = F_KK[p.first] * sum_KK;
        double N_Bp = F_Bp[p.first] * sum_Bp;
        double N_Bm = F_Bm[p.first] * sum_Bm;
        double N_Kspi0 = F_Kspi0[p.first] * sum_Kspi0;
        double N_K = K[p.first] * sumK;
        double N_Kbar = Kbar[p.first] * sumKbar;
        double Y_KK = Y_CP(p.first, 1, K0, Kbar0) * sum_KK;
        double YBp = Y_Bp(p.first, K0, Kbar0, Z0) * sum_Bp;
        double YBm = Y_Bm(p.first, K0, Kbar0, Z0) * sum_Bm;
        double Y_Kspi0 = Y_CP(p.first, -1, K0, Kbar0) * sum_Kspi0;
        INFO(p.first<<" "<<N_KK<<" "<<N_Kspi0<<" "<<N_Bp<<" "<<N_Bm);
        INFO(p.first<<" "<<Y_KK<<" "<<Y_Kspi0<<" "<<YBp<<" "<<YBm);
    }


    auto chi2_CP = [&F_KK, &F_Kspi0, sum_KK, sum_Kspi0, MPS, &Y_CP, &K, &Kbar](){
        double chi2 = 0;
        for (auto p : (*(&F_KK))){
            int bin = p.first;
            double N_KK = (*(&F_KK))[bin] * sum_KK;
            double N_Kspi0 = (*(&F_Kspi0))[bin] * sum_Kspi0;
            double E_KK = (*(&Y_CP))(bin, 1, (*(&K)), (*(&Kbar))) * sum_KK;
            double E_Kspi0 = (*(&Y_CP))(bin, -1, (*(&K)), (*(&Kbar))) * sum_Kspi0;
            double dN_KK = std::pow(N_KK, 0.5);
            double dE_KK = std::pow(E_KK, 0.5);
            double err_KK = std::pow(N_KK + E_KK, 0.5);

            double dN_Kspi0 = std::pow(N_Kspi0, 0.5);
            double dE_Kspi0 = std::pow(E_Kspi0, 0.5);
            double err_Kspi0 = std::pow(N_Kspi0 + E_Kspi0, 0.5);
            if (N_KK!=0)chi2 += std::pow((N_KK - E_KK)/dN_KK, 2);
            //chi2 += std::pow(N_KK - E_KK, 2)/err_KK;
           // pdf = mu^n/n! e^-mu
           // log pdf = n log mu - log n! - mu
            //chi2 += -2*(-E_KK + N_KK * std::log(E_KK) - std::log(std::tgamma(N_KK+1)));
            if(N_Kspi0!=0) chi2 += std::pow((N_Kspi0 - E_Kspi0)/dN_Kspi0, 2);
            //chi2 += std::pow(N_Kspi0 - E_Kspi0, 2)/err_Kspi0; 
            //chi2 +=-2*(-E_Kspi0 + N_Kspi0 * std::log(E_Kspi0) - std::log(std::tgamma(N_Kspi0+1)));
        }
        return chi2;
    };

    auto chi2_Kspipi = [&F_Kspipi, sum_Kspipi, &Y_Kspipi, &K, &Kbar, MPS](){
        double chi2 = 0;
        for (auto p : (*(&F_Kspipi))){
            std::pair<int, int> bin_pair = p.first;
            int bin_sig = bin_pair.first;
            int bin_tag = bin_pair.second;
            double N_Kspipi = p.second * sum_Kspipi;
            double dN_Kspipi = std::pow(N_Kspipi, 0.5);
            double E_Kspipi = (*(&Y_Kspipi))(bin_sig, bin_tag,  (*(&K)), (*(&Kbar))) * sum_Kspipi;
            double dE_Kspipi = std::pow(E_Kspipi, 0.5);
            double err_Kspipi = std::pow(N_Kspipi + E_Kspipi, 0.5);
            if (N_Kspipi!=0) chi2 += std::pow((E_Kspipi - N_Kspipi)/dN_Kspipi, 2);
             //chi2 += std::pow(E_Kspipi - N_Kspipi, 2)/err_Kspipi;

            //chi2 +=-2*( -E_Kspipi + N_Kspipi * std::log(E_Kspipi) - std::log(std::tgamma(N_Kspipi+1)));
            //if (N_Kspipi!=0) chi2 += 
        }
        return chi2;
    };

    auto chi2_BESIII = [&chi2_Kspipi, &chi2_CP, MPS](){
        double chi2 = 0;
        return (*(&chi2_Kspipi))() + (*(&chi2_CP))();
    };



    Minimiser mini_BESIII(chi2_BESIII, &MPS);
    mini_BESIII.doFit();
    for (int i=1;i<nBins + 1;i++){
        double ci_Fit = MPS["c_" + std::to_string(i)]->mean();
        double dci_Fit = MPS["c_" + std::to_string(i)]->err();
        double si_Fit = MPS["s_" + std::to_string(i)]->mean();
        double dsi_Fit = MPS["s_" + std::to_string(i)]->err();

        complex_t zi = Z0[i];
        double diff_ci = ci_Fit - std::real(zi);
        double diff_si = si_Fit - std::imag(zi);
        double pull_ci = diff_ci;
        double pull_si = diff_si;
        if (dci_Fit!= 0) pull_ci = pull_ci/dci_Fit;
        if (dsi_Fit!= 0) pull_si = pull_si/dsi_Fit;
        INFO(i<<" "<<std::real(zi)<<" "<<ci_Fit<<" +- "<<dci_Fit<<" "<<diff_ci<<" "<<pull_ci);
        INFO(i<<" "<<std::imag(zi)<<" "<<si_Fit<<" +- "<<dsi_Fit<<" "<<diff_si<<" "<<pull_si);
        

    }   

    auto chi2_B = [&F_Bp, &F_Bm, sum_Bp, sum_Bm, MPS, &Y_Bp ,&Y_Bm, &K0, &Kbar0, &Z0, nBins](){
        double chi2 = 0;
        std::map<int, complex_t> Z;
        for (int i =1;i<nBins + 1;i++){
            double c = MPS["c_"+std::to_string(i)]->mean();
            double s = MPS["s_"+std::to_string(i)]->mean();
            Z.insert(std::pair<int, complex_t>({i, complex_t(c, s)}));
            Z.insert(std::pair<int, complex_t>({-i, complex_t(c, -s)}));
        }
        for (auto p : (*(&F_Bp))){

            double _Y_Bp = MPS["YBp_"+std::to_string(p.first)]->mean();
            double _Y_Bm = MPS["YBm_"+std::to_string(p.first)]->mean();
            int bin = p.first;
            double N_Bp = (*(&F_Bp))[bin] * sum_Bp;
            double N_Bm = (*(&F_Bm))[bin] * sum_Bm;
            double E_Bp = (*(&Y_Bp))(bin, (*(&K0)), (*(&Kbar0)),(*(&Z))  ) * sum_Bp;
            double E_Bm = (*(&Y_Bm))(bin, (*(&K0)), (*(&Kbar0)),(*(&Z)) ) * sum_Bm;
            double dN_Bp = std::pow(N_Bp, 0.5);
            double dE_Bp = std::pow(E_Bp, 0.5);
            double err_Bp = std::pow(N_Bp + E_Bp, 0.5);

            double dN_Bm = std::pow(N_Bm, 0.5);
            double dE_Bm = std::pow(E_Bm, 0.5);
            double err_Bm = std::pow(N_Bm + E_Bm, 0.5);
            //if (N_Bp!=0)chi2 += std::pow((N_Bp - E_Bp), 2);
           if (N_Bp!=0)chi2 += std::pow((N_Bp - E_Bp)/dN_Bp, 2);
            //chi2 += std::pow(N_Bp - E_Bp, 2)/err_Bp;
           // pdf = mu^n/n! e^-mu
           // log pdf = n log mu - log n! - mu
            //chi2 +=  2 * E_Bp - 2 * N_Bp * std::log(E_Bp);// + 2 * std::log(std::tgamma(N_Bp + 1));
            //chi2 +=  2 * E_Bm - 2 * N_Bm * std::log(E_Bm);// + 2 * std::log(std::tgamma(N_Bm + 1));
            //if(N_Bm!=0) chi2 += std::pow((N_Bm - E_Bm)/dN_Bm, 2);
            if(N_Bm!=0) chi2 += std::pow((N_Bm - E_Bm)/dN_Bm, 2);
            //chi2 += std::pow(N_Bm - E_Bm, 2)/err_Bm; 
            //chi2 +=-2*(-E_Bm + N_Bm * std::log(E_Bm) - std::log(std::tgamma(N_Bm+1)));
        }
        return chi2;
    };
    INFO("chi2 = "<<chi2_B());
    Minimiser mini_B(chi2_B, &MPS);
    mini_B.doFit();
    double my_chi2 = 0;
    for (auto p : F_KK){
        double N_KK = F_KK[p.first] * sum_KK;
        double N_Bp = F_Bp[p.first] * sum_Bp;
        double N_Bm = F_Bm[p.first] * sum_Bm;
        double N_Kspi0 = F_Kspi0[p.first] * sum_Kspi0;
        double N_K = K[p.first] * sumK;
        double N_Kbar = Kbar[p.first] * sumKbar;
        double Y_KK = Y_CP(p.first, 1, K0, Kbar0) * sum_KK;
        double YBp = Y_Bp(p.first, K0, Kbar0, Z0) * sum_Bp;
        double YBm = Y_Bm(p.first, K0, Kbar0, Z0) * sum_Bm;
        double Y_Kspi0 = Y_CP(p.first, -1, K0, Kbar0) * sum_Kspi0;
        double pull_Bp = (N_Bp - YBp)/std::pow(N_Bp, 0.5);
        double pull_Bm = (N_Bm - YBm)/std::pow(N_Bm, 0.5);
        my_chi2 += std::pow(pull_Bp, 2) + std::pow(pull_Bm, 2);
        INFO(p.first<<" "<<N_KK<<" "<<N_Kspi0<<" "<<N_Bp<<" "<<N_Bm);
        INFO(p.first<<" "<<Y_KK<<" "<<Y_Kspi0<<" "<<YBp<<" "<<YBm<<" "<<pull_Bp<<" "<<pull_Bm);
        
    }
    INFO(chi2_B()<<" "<<my_chi2);

/*
    TFile * fScan = TFile::Open("MILLScan.root", "RECREATE");
    double xp0 = MPS["pCoherentSum::x+"]->mean();
    double yp0 = MPS["pCoherentSum::y+"]->mean();
    double xm0 = MPS["pCoherentSum::x-"]->mean();
    double ym0 = MPS["pCoherentSum::y-"]->mean();
    TGraph * g_xp = mini_B.scan(MPS["pCoherentSum::x+"], MPS["pCoherentSum::x+"]->mean() - 3 * MPS["pCoherentSum::x+"]->err(),MPS["pCoherentSum::x+"]->mean() + 3 * MPS["pCoherentSum::x+"]->err(), 0.1 * MPS["pCoherentSum::x+"]->err());

    MPS["pCoherentSum::x+"]->setCurrentFitVal(xp0);
    TGraph * g_yp = mini_B.scan(MPS["pCoherentSum::y+"], MPS["pCoherentSum::y+"]->mean() - 3 * MPS["pCoherentSum::y+"]->err(),MPS["pCoherentSum::y+"]->mean() + 3 * MPS["pCoherentSum::y+"]->err(), 0.1 * MPS["pCoherentSum::y+"]->err());   
    MPS["pCoherentSum::y+"]->setCurrentFitVal(yp0);

    TGraph * g_xm = mini_B.scan(MPS["pCoherentSum::x-"], MPS["pCoherentSum::x-"]->mean() - 3 * MPS["pCoherentSum::x-"]->err(),MPS["pCoherentSum::x-"]->mean() + 3 * MPS["pCoherentSum::x-"]->err(), 0.1 * MPS["pCoherentSum::x-"]->err());
    MPS["pCoherentSum::x-"]->setCurrentFitVal(xm0);
    TGraph * g_ym = mini_B.scan(MPS["pCoherentSum::y-"], MPS["pCoherentSum::y-"]->mean() - 3 * MPS["pCoherentSum::y-"]->err(),MPS["pCoherentSum::y-"]->mean() + 3 * MPS["pCoherentSum::y-"]->err(), 0.1 * MPS["pCoherentSum::y-"]->err());
    MPS["pCoherentSum::y-"]->setCurrentFitVal(ym0);




    g_xp->Write("xp");
    g_yp->Write("yp");
    g_ym->Write("ym");
    g_xm->Write("xm");
    fScan->Close();
   
    std::ofstream fYScan;
    fYScan.open("MILLplus_Scan.txt" ,std::ofstream::out | std::ofstream::app );

    int N = 10000;
    double x0 = -1;
    double xn = 1;
    int n =0 ;
    double dx = (xn - x0)/std::pow(N, 0.5);




    MPS["pCoherentSum::x+"]->fix();
    MPS["pCoherentSum::y+"]->fix();
    for (int i=0;i<(int)std::pow(N, 0.5);i++){
        for (int j=0;j<(int)std::pow(N, 0.5);j++){
            double _x = x0 + dx * i;
            double _y = x0 + dx * j;

            MPS["pCoherentSum::x+"]->setCurrentFitVal(_x);
            MPS["pCoherentSum::y+"]->setCurrentFitVal(_y);
            mini_B.doFit();
//            MPS["pCoherentSum::x-"]->setCurrentFitVal(_x);
//            MPS["pCoherentSum::y-"]->setCurrentFitVal(_y);

            fYScan<<_x<<" "<<_y<<" "<<chi2_B()<<"\n";
            
           
        }
    }
    fYScan.close();
    */
    

}