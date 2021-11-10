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
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");

  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;
auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 


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
  EventList mc =  Generator<>(eventType, &rndm).generate(NamedParameter<size_t>("NInt", 1e6));
  

  INFO("Testing a method to calc ACeif quickly using logs + exps");
  EventType KK({"D0", "K+", "K-"});
EventList mcKK =  Generator<>(KK, &rndm).generate(NamedParameter<size_t>("NInt", 1e6));
  pCorrelatedSum psi(eventType, KK, MPS);
  PhaseCorrection pc(MPS);


  psi.setEvents(mc, mcKK);
  psi.setMC(mc, mcKK);
  psi.prepare();
  ProfileClock clockManual;
  clockManual.start();
  double normManual = psi.norm();
  clockManual.stop();
  INFO("My value for norm by sum(|psi|^2) = "<<normManual<<" took "<<clockManual);
  return 0;

  CoherentSum A(eventType, MPS);
  CoherentSum B(KK, MPS);
  CoherentSum C(eventType.conj(true), MPS);
  CoherentSum D(KK.conj(true), MPS);
  A.setEvents(mc); A.setMC(mc); A.prepare();
  B.setEvents(mcKK); B.setMC(mcKK); B.prepare();
  C.setEvents(mc); C.setMC(mc); C.prepare();
  D.setEvents(mcKK); D.setMC(mcKK); D.prepare();

  std::vector<complex_t> As, Bs, Cs, Ds;
  for (int i=0;i<mc.size();i++){
    As.push_back(A.getVal(mc[i]));
    Bs.push_back(B.getVal(mcKK[i]));
    Cs.push_back(C.getVal(mc[i]));
    Ds.push_back(D.getVal(mcKK[i]));
  }

  ProfileClock clockcalcCorrL, clockcalcCorrXY;
  double fcalcCorrL = 0;
  double fcalcCorrXY = 0;
  clockcalcCorrL.start();

  #pragma omp parallel for reduction (+:fcalcCorrL)
  for (int i=0;i<mc.size();i++){
    double corr = pc.calcCorrL(mc[i]);
   fcalcCorrL +=std::norm(A.getVal(mc[i]) * B.getVal(mcKK[i])  * exp(complex_t(0, corr/2))- C.getVal(mc[i])  * D.getVal(mcKK[i]) * exp(complex_t(0, -corr/2)));
   //fcalcCorrL +=std::norm(A.getVal(mc[i]) * B.getVal(mcKK[i]) - C.getVal(mc[i])  * D.getVal(mcKK[i]) );


  }
  clockcalcCorrL.stop();
  std::vector<double> s01(mc.size()), s02(mc.size());

  for (int i=0;i<mc.size();i++){
      s01[i] = mc[i].s(0,1);
      s02[i] = mc[i].s(0,2);
  }

  

  clockcalcCorrXY.start();
  std::vector<double> corrXY(mc.size());
  #pragma omp parallel for
  for (int i=0;i<mc.size();i++){
    double _corrXY = pc.calcCorrXY(s01[i], s02[i]);
    corrXY[i] = std::norm(As[i] * Bs[i] * exp(complex_t(0, _corrXY/2)) - Cs[i] * Ds[i] * exp(complex_t(0, -_corrXY/2)));
  }


  for (int i=0;i<mc.size();i++){

    //fcalcCorrXY += std::norm(As[i] * Bs[i] * exp(complex_t(0, corrXY[i]/2)) - Cs[i] * Ds[i] * exp(complex_t(0, -corrXY[i]/2)));
    //fcalcCorrXY += std::norm(As[i] * Bs[i] - Cs[i] * Ds[i]);
    fcalcCorrXY += corrXY[i];

  }
  clockcalcCorrXY.stop();

  INFO("Got "<<fcalcCorrL<<" took "<<clockcalcCorrL<<" for f(evt)");
  INFO("Got "<<fcalcCorrXY<<" took "<<clockcalcCorrXY<<" for f(s01, s02)");



  ProfileClock clockgetVal;
  ProfileClock clockgetValForNorm;
  complex_t v0=0;
  complex_t v0forNorm=0;
  std::vector<double> timesGetVal;
  std::vector<double> timesGetValForNorm;
  for (int i=0;i<mc.size();i++){
  clockgetVal.start();
  v0 = psi.getVal(mc[i], mcKK[i]);
  clockgetVal.stop();
  timesGetVal.push_back(clockgetVal);
  }



  for (int i=0;i<mc.size();i++){
  clockgetValForNorm.start();
  v0forNorm = psi.getValForNorm(i);
  clockgetValForNorm.stop();
  
  timesGetValForNorm.push_back(clockgetValForNorm);
  }

  double meanGetValTime = std::accumulate(timesGetVal.begin(), timesGetVal.end(), 0);
  double meanGetValForNormTime = std::accumulate(timesGetValForNorm.begin(), timesGetValForNorm.end(), 0);
  INFO("mean t(getVal) = "<<meanGetValTime);
  INFO("mean t(getValForNorm) = "<<meanGetValForNormTime);

  INFO("v0 = "<<v0<<" took "<<clockgetVal);
  INFO("v0 = "<<v0forNorm<<" took "<<clockgetValForNorm);


  int nL=0;
  for (int i=0;i<mc.size();i++){
clockgetValForNorm.start();
  v0forNorm = psi.getValForNorm(i);
  clockgetValForNorm.stop();
   clockgetVal.start();
  v0 = psi.getVal(mc[i], mcKK[i]);
  clockgetVal.stop();
  if (clockgetValForNorm > clockgetVal){
    nL+=1;
  }
  
  }
INFO("Why did "<<nL<<" take longer?");
 


  INFO("log10 = "<<log(10));
  double x0 = 3;
  INFO("leg_8(x) = "<<std::legendre(8, x0));

  real_t N = mc.size();
  real_t nA=0;
  real_t nC = 0;
  ProfileClock clockNormSum;
  clockNormSum.start();
  double normSum = 0;
  for (int i=0;i<mc.size();i++){
    normSum += std::norm(psi.getVal(mc[i], mcKK[i]))/N;
  }
  clockNormSum.stop();
  INFO("normSum = "<<normSum<<" took "<<clockNormSum);
  std::vector<double> psi2Arr;
  ProfileClock clockPsi2Arr;
  ProfileClock clockCorr;


  std::vector<double> absA, absC, cosdd, sindd, dd;
  std::map<std::pair<size_t, size_t>, std::vector<double> > Pmn;
  std::map<std::pair<size_t, size_t>, MinuitParameter*> Cmn;
  
  double avg_dd = 0;
  ProfileClock clockAmp;
  clockAmp.start();
  for (auto evt:mc){
    complex_t a = A.getValNoCache(evt);
    complex_t c = C.getValNoCache(evt);
    absA.push_back(std::abs(a));
    absC.push_back(std::abs(c));
    nA += std::norm(a);
    nC += std::norm(c);
    double _dd = std::arg(a * std::conj(c));
    dd.push_back(_dd);
    cosdd.push_back(cos(_dd));
    sindd.push_back(sin(_dd));
    avg_dd += _dd/N;
    for (auto p : Pmn){
       double poly2d = pc.Poly2D(pc.getXY(evt)[0], pc.getXY(evt)[1], p.first.first, p.first.second);
       Pmn[p.first].push_back(poly2d);
    }
  }
  clockAmp.stop();
  INFO("average dd = "<<avg_dd<<" took "<<clockAmp);
  INFO("Cmn has "<<Cmn.size()<<" elements, Pmn has "<<Pmn.size()<<" elements ");

  ProfileClock clockMyNorm;
  clockMyNorm.start();
  clockPsi2Arr.start();
  for (int i=0;i<mc.size();i++){
  clockCorr.start();
    double corr = pc.calcCorrL(mc[i]);
  clockCorr.stop();
    //double sum = absA[i] * absC[i] * cos(dd[i] + corr);
    double sum = std::real(A.getVal(mc[i]) * std::conj(C.getVal(mc[i])) * exp(complex_t(0, corr))); 

    psi2Arr.push_back(sum);
  } 
  clockPsi2Arr.stop();
  INFO("Took "<<clockPsi2Arr<<" to vectorise the psis");
  INFO("Took "<<clockCorr<<" to calculate correction");

  ProfileClock clockAcc;
  clockAcc.start();
  double myNorm = std::accumulate(psi2Arr.begin(), psi2Arr.end(),0);
  clockAcc.stop();

  myNorm = (nA + nC - 2 * myNorm)/N;
  clockMyNorm.stop();

  INFO("Took "<<clockAcc<<" to acculmate");
  INFO("Got "<<myNorm<<" for normalisation took "<<clockMyNorm);
  double myNormLoop = 0;
  ProfileClock clockLoop;
  clockLoop.start();
  for (int i=0;i<psi2Arr.size();i++){
    myNormLoop += psi2Arr[i];
  }
  clockLoop.stop();
  INFO("Took "<<clockLoop<<" to loop");


  ProfileClock clockgetC;
  clockgetC.start();
  double c01 = MPS["PhaseCorrection::C0_1"]->mean();
  clockgetC.stop();
  INFO("Took "<<clockgetC<<" to get 1 MPS param");

  return 0;
 



}