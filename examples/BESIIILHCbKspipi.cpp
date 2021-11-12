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
#include "AmpGen/Minimiser.h"
#include "AmpGen/ProfileClock.h"
#include <TMath.h>
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CombGamLL.h"
#include "AmpGen/ProgressBar.h"
#include "AmpGen/FitResult.h"
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
#include <cmath>
using namespace AmpGen;
using namespace std::complex_literals;


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


  size_t Order = NamedParameter<size_t>("PhaseCorrection::Order", 2);
  std::string BESIIIDataFile = NamedParameter<std::string>("BESIIIDataSample", ""          , "Name of file containing data sample to fit." );
  std::string LHCbDataFile = NamedParameter<std::string>("LHCbDataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");

  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  auto BESIIITags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();  

  auto LHCbTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();

  size_t NInt = NamedParameter<size_t>("NInt", 1e6);

  INFO("Combined Fit for BESIII + LHCb");




  if( BESIIIDataFile == "" || LHCbDataFile== "" ) FATAL("Must specify input with option " << italic_on << "BESIIIDataSample" << italic_off <<" and "<<italic_on<< " LHCbDataSample"<<italic_off);
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);
  if (intFile == ""){

  }
  TRandom3 rndmSig, rndmTag;
  rndmSig.SetSeed( seed + 100 );
  rndmTag.SetSeed( seed + 100 );


  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif


 MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
//    add_CP_conjugate(MPS);
      AddCPConjugate(MPS);
  }

  EventType sigType = EventType(pNames);
  gRandom = &rndmSig;
  EventList mcSig =  Generator<>(sigType, &rndmSig).generate(NInt);
  std::vector<complex_t> ac (NInt);
  PhaseCorrection pc(MPS);
  CoherentSum A (sigType, MPS);
  CoherentSum C (sigType.conj(true), MPS);
  
  A.setEvents(mcSig); A.setMC(mcSig); A.prepare();
  C.setEvents(mcSig); C.setMC(mcSig); C.prepare();
  real_t normA = A.norm();
  real_t normC = C.norm();
  for (size_t i=0;i<mcSig.size();i++){
    ac[i] = A.getValNoCache(mcSig[i]) * std::conj(C.getValNoCache(mcSig[i]));
  }
  INFO("Calculating the int AC eif part now and share it between amps!");
  auto ACeif = [&pc, &ac, &mcSig](){
    
    real_t x=0; 
    real_t y=0; 
    #pragma omp parallel for reduction(+:x,y)
    for (size_t i=0;i<mcSig.size();i++){
      complex_t r = exp(complex_t(0,pc.calcCorrL(mcSig[i]))) * ac[i];
      x += std::real(r);
      y += std::imag(r);
    }
    complex_t z(x, y);
    return z/(real_t)mcSig.size();
  };
  complex_t aceif = ACeif();
  std::map<std::string, EventList> LHCbSig;
  for (auto tag : LHCbTags){
    auto arr_tag = split(tag, ' ');
    auto tagName = arr_tag[0];

    EventList Events = EventList(LHCbDataFile + ":" + tagName, sigType);
    LHCbSig.insert(std::pair<std::string, EventList>({tagName, Events}));
  }
  std::vector<pCoherentSum> LHCbAmps(LHCbSig.size());
  int i=0;
  INFO("Preparing B Amplitudes");
  for (auto tag : LHCbTags){
    auto arr_tag = split(tag, ' ');
    auto tagName = arr_tag[0];
    bool conjD =(bool) std::stoi(arr_tag[2]);
    int BSign = std::stoi(arr_tag[3]);
    bool useXY = (bool) std::stoi(arr_tag[4]);

    LHCbAmps[i] = pCoherentSum(sigType, MPS, BSign, useXY, conjD);
    LHCbAmps[i].setEventsByRef(&LHCbSig[tagName]);
    LHCbAmps[i].setMCByRef(&mcSig);
    LHCbAmps[i].updateZACst(aceif); 
    INFO("Preparing LHCb "<<tagName);
    ProfileClock clockPreparePsiLHCb;
    clockPreparePsiLHCb.start();
    LHCbAmps[i].prepare();
    clockPreparePsiLHCb.stop();
    
    INFO("Prepared "<<tagName<<" amp in "<<clockPreparePsiLHCb); 


    ProfileClock clockNormPsiLHCb;
    clockNormPsiLHCb.start();
    double normPsi = LHCbAmps[i].norm();
    clockNormPsiLHCb.stop();
    INFO("Norm for "<<tagName <<" = "<<normPsi<<" took "<<clockNormPsiLHCb);
    
    i++;
    
  }

  //std::map<std::string, pCorrelatedSum> BESIIIAmps;
  std::map<std::string, EventList> BESIIISig;
  std::map<std::string, EventList> BESIIITag;
  std::map<std::string, EventList> BESIIITagMC;

  std::vector< complex_t> BDs;

  INFO("Signal events = "<<sigType<<" have "<<mcSig.size()<<" MC events");
  for(auto tag : BESIIITags){
    auto arr_tag = split(tag, ' ');
    auto tagName = arr_tag[0];
    EventType tagType = Particle(arr_tag[1], {}, false).eventType();
    EventList sigEvents = EventList(BESIIIDataFile + ":Signal_" + tagName, sigType);
    EventList tagEvents = EventList(BESIIIDataFile + ":Tag_" + tagName, tagType); 
    INFO("Tag Type = "<<tagType);
    
  gRandom = &rndmSig;
    EventList mcTag = Generator<>(tagType, &rndmSig).generate(NInt);
    BESIIITagMC.insert(std::pair<std::string, EventList>(tagName, mcTag));
    BESIIISig.insert(std::pair<std::string, EventList>(tagName, sigEvents));
    BESIIITag.insert(std::pair<std::string, EventList>(tagName, tagEvents));
  }

  //Hard code our tags for now 

  std::vector<pCorrelatedSum> BESIIIAmps (BESIIITags.size());
  i =0 ;
  for (auto p:BESIIISig){
    EventList sig = BESIIISig[p.first];
    EventList tag = BESIIITag[p.first];
    EventList tagMC = BESIIITagMC[p.first];
    BESIIIAmps[i] = pCorrelatedSum (sig.eventType(), tag.eventType(), MPS);
    BESIIIAmps[i].setEventsByRef(&sig, &tag);
  
    BESIIIAmps[i].setMCByRef(&mcSig, &tagMC);
    BESIIIAmps[i].updateZACst(aceif);
    
    BESIIIAmps[i].prepare();
    CoherentSum B(tag.eventType(), MPS); B.setEvents(tagMC); B.setMC(tagMC); B.prepare();
    CoherentSum D(tag.eventType().conj(true), MPS); D.setEvents(tagMC); D.setMC(tagMC); D.prepare();
    real_t normB = B.norm();
    real_t normD = D.norm();

    real_t norm = BESIIIAmps[i].norm();
    INFO("norm = "<<norm);
    INFO("normA = "<<normA);
    INFO("normB = "<<normB);
    INFO("normC = "<<normC);
    INFO("normD = "<<normD);

    INFO("Calculating norm manually too");
    real_t normMan = 0;

    INFO("sigMC s(0,1)[0] = "<<mcSig[0].s(0,1));
    INFO("tagMC s(0,1)[0] = "<<tagMC[0].s(0,1));
    complex_t BD = BESIIIAmps[i].getBDstSum();
    BDs.push_back(BD);
 real_t myNormCalc = BESIIIAmps[i].normFromZ(aceif, BD);
INFO("my attempt at getting norm  = "<<myNormCalc);
    INFO("BD = "<<BD<<" ACeif = "<<aceif);
    complex_t zAC = 0;
    real_t N = mcSig.size();

 

    
//    ps[i].debugNorm();
    i++;
  }
  


  //  complex_t zAC = 0;
    real_t N = mcSig.size();


  ProfileClock clockAC, clockNorm;

  clockNorm.start();
  clockAC.start();
  /*
    for (size_t j=0;j<mcSig.size();j++){
      real_t corr = pc.calcCorrL(mcSig[j]);   
      zAC += A.getValNoCache(mcSig[j]) * std::conj(C.getValNoCache(mcSig[j])) * exp(complex_t(0, corr))/N;
    }
    */
   complex_t zAC = ACeif();
  clockAC.stop();

  ProfileClock clockBESIII;
  clockBESIII.start();
  for (auto p : BESIIIAmps){
    p.updateZACst(zAC);
    INFO("n = "<<p.norm());
  }
  clockBESIII.stop();
  




  ProfileClock clockLHCb;
  clockLHCb.start();

  for (auto p : LHCbAmps){
    p.updateZACst(zAC);
    INFO("n = "<<p.norm());
  }

  clockLHCb.stop();

clockNorm.stop();
  INFO("Took "<<clockAC<<" to calc int ACeif");
  INFO("Took "<<clockBESIII<<" to calc BESIII Norms");
  INFO("Took "<<clockLHCb<<" to calc LHCb Norms");
  INFO("Took "<<clockNorm<<" to calculate all "<<BESIIIAmps.size()<<" normalisations for BESIII + "<<LHCbAmps.size()<<" normalisations for LHCb");

  auto LL_BESIII = [&BESIIIAmps, &ACeif, &BESIIISig, &BESIIITag](){
    complex_t z = ACeif();
    real_t ll = 0;
    int i=0;
    for (auto p : BESIIISig){
      auto evts = p.second;
      BESIIIAmps[i].updateZACst(z);
      real_t ll_i=0;
      real_t n = BESIIIAmps[i].norm();
      #pragma omp parallel for reduction (+:ll_i)
      for(size_t j=0;j<evts.size();j++){
        ll_i += log(std::norm(BESIIIAmps[i].getValNoCache(BESIIISig[p.first][j], BESIIITag[p.first][j]))/n);
      }
      ll+=ll_i;
      i++;
//      ll += p.LL();

    }
    return -2*ll;
  };

  auto LL_LHCb = [&LHCbAmps, &ACeif, &LHCbSig](){
    complex_t z = ACeif();
    real_t ll = 0;
    int i=0;
    for (auto p : LHCbSig){
      auto evts = p.second;
      LHCbAmps[i].updateZACst(z);
      real_t ll_i=0;
      real_t n = LHCbAmps[i].norm();
      #pragma omp parallel for reduction (+:ll_i)
      for(size_t j=0;j<evts.size();j++){
        ll_i += log(std::norm(LHCbAmps[i].getValNoCache(LHCbSig[p.first][j]))/n);
      }
      ll+=ll_i;
      i++;
//      ll += p.LL();

    }
    return -2*ll;
  };
  auto LL_BESIII_LHCb = [&BESIIIAmps, &LHCbAmps, &ACeif, &BESIIISig, &BESIIITag, &LHCbSig](){
    complex_t z = ACeif();
    real_t ll = 0;
    int i=0;
    for (auto p : BESIIISig){
      auto evts = p.second;
      BESIIIAmps[i].updateZACst(z);
      real_t ll_i=0;
      real_t n = BESIIIAmps[i].norm();
      #pragma omp parallel for reduction (+:ll_i)
      for(size_t j=0;j<evts.size();j++){
        ll_i += log(std::norm(BESIIIAmps[i].getValNoCache(BESIIISig[p.first][j], BESIIITag[p.first][j]))/n);
      }
      ll+=ll_i;
      i++;
//      ll += p.LL();

    }
    i = 0;
    for (auto p : LHCbSig){
      auto evts = p.second;
      LHCbAmps[i].updateZACst(z);
      real_t ll_i=0;
      real_t n = LHCbAmps[i].norm();
      #pragma omp parallel for reduction (+:ll_i)
      for(size_t j=0;j<evts.size();j++){
        ll_i += log(std::norm(LHCbAmps[i].getValNoCache(LHCbSig[p.first][j]))/n);
      }
      ll+=ll_i;
      i++;
//      ll += p.LL();

    }
    return -2*ll;
  };

  auto LL_BESIII_LHCb_LASSO = [&LL_BESIII_LHCb, &MPS, &Order](){

        real_t lambda = MPS["LASSO::lambda"]->mean() ;
        real_t penalty=0;
        for (size_t i=0;i<Order+1;i++){
            for (size_t j=0;j<Order+1-i;j++){
                int i1 = i;
                int i2 = 2 * j + 1;
                auto p = MPS["PhaseCorrection::C" + std::to_string(i1) + "_" + std::to_string(i2)];
                penalty += std::abs(p->mean());
            }
        }
        return LL_BESIII_LHCb() + lambda * penalty;
  };





    ProfileClock clockBESIIILL;
    clockBESIIILL.start();
    INFO("LL = "<<LL_BESIII());
    clockBESIIILL.stop();
    INFO("Took "<<clockBESIIILL<<" to calc LL(BESIII)");

    ProfileClock clockLHCbLL;
    clockLHCbLL.start();
    INFO("LL = "<<LL_LHCb());
    clockLHCbLL.stop();
    INFO("Took "<<clockLHCbLL<<" to calc LL(LHCb)");

    ProfileClock clockBESIII_LHCbLL;
    clockBESIII_LHCbLL.start();
    INFO("LL = "<<LL_BESIII_LHCb());
    clockBESIII_LHCbLL.stop();
    INFO("Took "<<clockBESIII_LHCbLL<<" to calc LL(BESIII_LHCb)");

    
    size_t maxCalls       = NamedParameter<size_t>( "Minimiser::MaxCalls"  , 100000);
    real_t maxTime_ms = maxCalls * clockBESIII_LHCbLL;
    INFO("Max time a fit should take = "<<maxTime_ms<<" ms = "<<maxTime_ms/(60 * 60 * 1e3)<<" hours");



    real_t lambda = MPS["LASSO::lambda"]->mean() ;
    if (lambda != 0){
      INFO("Doing LASSO Fit with lambda = "<<lambda);
      Minimiser mini_BESIII_LHCb_LASSO(LL_BESIII_LHCb_LASSO, &MPS);
      mini_BESIII_LHCb_LASSO.doFit();
      real_t LASSO_Threshold = NamedParameter<real_t>("LASSO::Threshold", 0.01);
      for (long unsigned int i=0;i<Order+1;i++){
        for (long unsigned int j=0;j<Order+1-i;j++){
            auto p = MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j + 1)];
            if (abs(p->mean())< LASSO_Threshold){
              MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j + 1)]->setCurrentFitVal(0);
              MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j + 1)]->fix();
              INFO("Excluding C"<<i<<"_"<<2*j+1<<" from fit");
            }

        }
      }
    }
    Minimiser mini_BESIII_LHCb(LL_BESIII_LHCb, &MPS);
    mini_BESIII_LHCb.doFit();
    FitResult * fr = new FitResult(mini_BESIII_LHCb);
    fr->writeToFile(logFile);

        








 
  return 0;
  INFO("Out of loop i = "<<i);






  /*
  std::map<std::string, pCoherentSum> LHCbAmps;
  for (auto tag : LHCbTags){
    auto arr_tag = split(tag, ' ');
    auto tagName = arr_tag[0];
    bool conjD =(bool) std::stoi(arr_tag[2]);
    int BSign = std::stoi(arr_tag[3]);
    bool useXY = (bool) std::stoi(arr_tag[4]);
    EventList Events = EventList(LHCbDataFile + ":" + tagName, sigType);
    pCoherentSum psi(sigType, MPS, BSign, useXY, conjD);
    psi.setEvents(Events);
    psi.setMC(mcSig);
    ProfileClock clockPreparePsiLHCb;
    clockPreparePsiLHCb.start();
    psi.prepare();
    clockPreparePsiLHCb.stop();
    
    INFO("Prepared "<<tagName<<" amp in "<<clockPreparePsiLHCb); 

    LHCbAmps.insert(std::pair<std::string, pCoherentSum>({tagName, psi}));
    ProfileClock clockNormPsiLHCb;
    clockNormPsiLHCb.start();
    double normPsi = psi.norm();
    clockNormPsiLHCb.stop();
    INFO("Norm for "<<tagName <<" = "<<normPsi<<" took "<<clockNormPsiLHCb);

    
  }

  auto LL_BESIII = [&BESIIIAmps](){
    double LL = 0;
    for (auto p : BESIIIAmps){
      LL += p.second.LL();
    }
    return LL;
  };

  auto LL_LHCb = [&LHCbAmps](){
    double LL =0;
    for(auto p : LHCbAmps){
      LL += p.second.LL();
    }
    return LL;
  };

  ProfileClock clockLL_BESIII, clockLL_LHCb;
  clockLL_BESIII.start();
  complex_t LL_BESIII_0 = BESIIIAmps["Kspipi"].getVal(mcSig[0], mcSig[0]);
  clockLL_BESIII.stop();
  INFO("Got LL_BESIII = "<<LL_BESIII_0<<" in "<<clockLL_BESIII);


  clockLL_LHCb.start();
  complex_t LL_LHCb_0 = LHCbAmps["Bm2Dhm"].getVal(mcSig[0]);
  clockLL_LHCb.stop();

  INFO("Got LL_LHCb = "<<LL_LHCb_0<<" in "<<clockLL_LHCb);


*/


  return 0;
}