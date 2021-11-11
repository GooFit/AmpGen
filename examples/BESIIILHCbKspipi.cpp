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
    complex_t z = 0;
    for (size_t i=0;i<mcSig.size();i++){
      z += exp(complex_t(0,pc.calcCorrL(mcSig[i]))) * ac[i];
    }
    return z/(real_t)mcSig.size();
  };
  complex_t aceif = ACeif();
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

  std::vector<pCorrelatedSum> ps (BESIIITags.size());
  size_t i =0 ;
  for (auto p:BESIIISig){
    EventList sig = BESIIISig[p.first];
    EventList tag = BESIIITag[p.first];
    EventList tagMC = BESIIITagMC[p.first];
    ps[i] = pCorrelatedSum (sig.eventType(), tag.eventType(), MPS);
    ps[i].setEventsByRef(&sig, &tag);
  
    ps[i].setMCByRef(&mcSig, &tagMC);

    
    ps[i].prepare();
    CoherentSum B(tag.eventType(), MPS); B.setEvents(tagMC); B.setMC(tagMC); B.prepare();
    CoherentSum D(tag.eventType().conj(true), MPS); D.setEvents(tagMC); D.setMC(tagMC); D.prepare();
    real_t normB = B.norm();
    real_t normD = D.norm();

    real_t norm = ps[i].norm();
    INFO("norm = "<<norm);
    INFO("normA = "<<normA);
    INFO("normB = "<<normB);
    INFO("normC = "<<normC);
    INFO("normD = "<<normD);

    INFO("Calculating norm manually too");
    real_t normMan = 0;

    INFO("sigMC s(0,1)[0] = "<<mcSig[0].s(0,1));
    INFO("tagMC s(0,1)[0] = "<<tagMC[0].s(0,1));
    complex_t BD = ps[i].getBDstSum();
    BDs.push_back(BD);
 real_t myNormCalc = ps[i].normFromZ(aceif, BD);
INFO("my attempt at getting norm  = "<<myNormCalc);
    INFO("BD = "<<BD<<" ACeif = "<<aceif);
    real_t nA = 0;
    real_t nB = 0;
    real_t nC = 0;
    real_t nD = 0;
    complex_t zAC = 0;
    complex_t zBD = 0;
    real_t nAB = 0;
    real_t nCD = 0;
    complex_t zABCD = 0;
    real_t N = mcSig.size();
    for (size_t j=0;j<mcSig.size();j++){
   
      nA += std::norm(A.getValNoCache(mcSig[j]))/N;
      nB += std::norm(B.getValNoCache(tagMC[j]))/N;
      nC += std::norm(C.getValNoCache(mcSig[j]))/N;
      nD += std::norm(D.getValNoCache(tagMC[j]))/N;
      zAC += A.getValNoCache(mcSig[j]) * std::conj(C.getValNoCache(mcSig[j]))/N;
      zBD += B.getValNoCache(tagMC[j]) * std::conj(D.getValNoCache(tagMC[j]))/N;

    }
    INFO("manualA = "<<nA);
    INFO("manualB = "<<nB);
    INFO("manualC = "<<nC);
    INFO("manualD = "<<nD);
    INFO("zAC = "<<zAC);
    INFO("zBD = "<<zBD);
    
    normMan = nA * nB + nC * nD - 2 * std::real(zAC * zBD);
    INFO("normManual = "<<normMan);
    //if (sig.eventType()==tag.eventType()) BD = std::conj(aceif);
   

    
//    ps[i].debugNorm();
    i++;
  }
  return 0;
  INFO("Out of loop i = "<<i);
  
  ProfileClock clockNormKK;
  ProfileClock clockLLKK;
  clockLLKK.start();
  clockNormKK.start();
  double psKKNorm = ps[0].norm();
  clockNormKK.stop();
  double LLKK = 0; 
  for (size_t i=0;i<BESIIISig["KK"].size();i++){
    Event sig = BESIIISig["KK"][i];
    Event tag = BESIIITag["KK"][i];
    LLKK += log(std::norm(ps[0].getValNoCache(sig, tag))/psKKNorm );
 
  }
  clockLLKK.stop();
  INFO("norm0 = "<<psKKNorm<<" took "<<clockNormKK);
  INFO("LLKK = "<<LLKK<<" took "<<clockLLKK);

  ProfileClock clockNormKspi0;
  ProfileClock clockLLKspi0;
  clockLLKspi0.start();
  clockNormKspi0.start();
  double psKspi0Norm = ps[1].norm();
  clockNormKspi0.stop();
  double LLKspi0 = 0; 
  for (size_t i=0;i<BESIIISig["Kspi0"].size();i++){
    Event sig = BESIIISig["Kspi0"][i];
    Event tag = BESIIITag["Kspi0"][i];
    LLKspi0 += log(std::norm(ps[1].getValNoCache(sig, tag))/psKspi0Norm );
 
  }
  clockLLKspi0.stop();
  INFO("norm0 = "<<psKspi0Norm<<" took "<<clockNormKspi0);
  INFO("LLKspi0 = "<<LLKspi0<<" took "<<clockLLKspi0);

  ProfileClock clockNormKppim;
  ProfileClock clockLLKppim;
  clockLLKppim.start();
  clockNormKppim.start();
  double psKppimNorm = ps[2].norm();
  clockNormKppim.stop();
  double LLKppim = 0; 
  for (size_t i=0;i<BESIIISig["Kppim"].size();i++){
    Event sig = BESIIISig["Kppim"][i];
    Event tag = BESIIITag["Kppim"][i];
    LLKppim += log(std::norm(ps[2].getValNoCache(sig, tag))/psKppimNorm );
 
  }
  clockLLKppim.stop();
  INFO("norm0 = "<<psKppimNorm<<" took "<<clockNormKppim);
  INFO("LLKppim = "<<LLKppim<<" took "<<clockLLKppim);

  ProfileClock clockNormKmpip;
  ProfileClock clockLLKmpip;
  clockLLKmpip.start();
  clockNormKmpip.start();
  double psKmpipNorm = ps[3].norm();
  clockNormKmpip.stop();
  double LLKmpip = 0; 
  for (size_t i=0;i<BESIIISig["Kmpip"].size();i++){
    Event sig = BESIIISig["Kmpip"][i];
    Event tag = BESIIITag["Kmpip"][i];
    LLKmpip += log(std::norm(ps[3].getValNoCache(sig, tag))/psKmpipNorm );
 
  }
  clockLLKmpip.stop();
  INFO("norm0 = "<<psKmpipNorm<<" took "<<clockNormKmpip);
  INFO("LLKmpip = "<<LLKmpip<<" took "<<clockLLKmpip);

  ProfileClock clockNormKspipi;
  ProfileClock clockLLKspipi;

  Event sig = BESIIISig["Kspipi"][0];
  Event tag = BESIIITag["Kspipi"][0];

  complex_t v_kspipi_00 = ps[4].getValNoCache(sig, tag);
  INFO("vKspipi "<<v_kspipi_00);


  clockLLKspipi.start();
  clockNormKspipi.start();

 
  double psKspipiNorm = ps[4].norm();
  clockNormKspipi.stop();
 return 0;
  double LLKspipi = 0; 
  for (size_t i=0;i<BESIIISig["Kspipi"].size();i++){
    Event sig = BESIIISig["Kspipi"][i];
    Event tag = BESIIITag["Kspipi"][i];
    LLKspipi += log(std::norm(ps[4].getValNoCache(sig, tag))/psKspipiNorm );
 
  }
  clockLLKspipi.stop();
  INFO("norm0 = "<<psKspipiNorm<<" took "<<clockNormKspipi);
  INFO("LLKspipi = "<<LLKspipi<<" took "<<clockLLKspipi);







  return 0;

  EventList sigKK = BESIIISig["KK"];
  EventList tagKK = BESIIITag["KK"];
  EventList tagKKMC = BESIIITagMC["KK"];
  pCorrelatedSum psKK(sigKK.eventType(), tagKK.eventType(), MPS);
  psKK.setEvents(sigKK, tagKK);
  psKK.setMC(mcSig, tagKKMC);
  psKK.prepare();
  real_t normKK = psKK.norm();
  INFO("normKK = "<<normKK);

  EventList sigKspi0 = BESIIISig["Kspi0"];
  EventList tagKspi0 = BESIIITag["Kspi0"];
  EventList tagKspi0MC = BESIIITagMC["Kspi0"];
  pCorrelatedSum psKspi0(sigKspi0.eventType(), tagKspi0.eventType(), MPS);
  psKspi0.setEvents(sigKspi0, tagKspi0);
  psKspi0.setMC(mcSig, tagKspi0MC);
  psKspi0.prepare();
  real_t normKspi0 = psKspi0.norm();
  INFO("normKspi0 = "<<normKspi0);

  EventList sigKmpip = BESIIISig["Kmpip"];
  EventList tagKmpip = BESIIITag["Kmpip"];
  EventList tagKmpipMC = BESIIITagMC["Kmpip"];
  pCorrelatedSum psKmpip(sigKmpip.eventType(), tagKmpip.eventType(), MPS);
  psKmpip.setEvents(sigKmpip, tagKmpip);
  psKmpip.setMC(mcSig, tagKmpipMC);
  psKmpip.prepare();
  real_t normKmpip = psKmpip.norm();
  INFO("normKmpip = "<<normKmpip);

  EventList sigKppim = BESIIISig["Kppim"];
  EventList tagKppim = BESIIITag["Kppim"];
  EventList tagKppimMC = BESIIITagMC["Kppim"];
  pCorrelatedSum psKppim(sigKppim.eventType(), tagKppim.eventType(), MPS);
  psKppim.setEvents(sigKppim, tagKppim);
  psKppim.setMC(mcSig, tagKppimMC);
  psKppim.prepare();
  real_t normKppim = psKppim.norm();
  INFO("normKppim = "<<normKppim);

  EventList sigKspipi = BESIIISig["Kspipi"];
  EventList tagKspipi = BESIIITag["Kspipi"];

  pCorrelatedSum psKspipi(sigKspipi.eventType(), tagKspipi.eventType(), MPS);
  psKspipi.setEvents(sigKspipi, tagKspipi);
  psKspipi.setMC(mcSig, mcSig);
  psKspipi.prepare();
  real_t normKspipi = psKspipi.norm();
  INFO("normKspipi = "<<normKspipi);

  auto LL_BESIII = [&psKK, &psKspi0](){
    return psKK.LL() + psKspi0.LL();
  };

  double LL = LL_BESIII();







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