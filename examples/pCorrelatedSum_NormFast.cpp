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

  CoherentSum A(eventType, MPS);
  CoherentSum C(eventType.conj(true), MPS);
  A.setEvents(mc); A.setMC(mc); A.prepare();
  C.setEvents(mc); C.setMC(mc); C.prepare();

  INFO("log10 = "<<log(10));
  ProfileClock clockManual;
  clockManual.start();
  double normManual = psi.norm();
  clockManual.stop();
  INFO("My value for norm by sum(|psi|^2) = "<<normManual<<" took "<<clockManual);

  double nA = 0;
  double nC = 0;
  ProfileClock clockAC;
  clockAC.start();
  for (size_t i=0;i<mc.size();i++){
      nA += std::norm(A.getValNoCache(mc[i]));
      nC += std::norm(A.getValNoCache(mc[i]));
  }
  nA = nA/(real_t)mc.size();
  nC = nC/(real_t)mc.size();
  clockAC.stop();
  INFO("nA = "<<nA<<" nC = "<<nC<<" took "<<clockAC);


  std::vector<complex_t> g1;
  complex_t G1 = 0;

  std::map<std::pair<size_t, size_t>, std::vector<double> > g2;
  std::map<std::pair<size_t, size_t>, double> G2;

  std::vector<complex_t> Delta1;
  std::map<std::pair<size_t, size_t>, std::vector<double> > Delta2;
  real_t N = mc.size(); 
  ProfileClock clockG1;
  clockG1.start();
  double sumDD = 0;
  for (size_t i=0;i<mc.size();i++){
      complex_t a = A.getVal(mc[i]);
      complex_t c = C.getVal(mc[i]);
      double dd = std::arg(a * c);
      sumDD += dd;
      complex_t _g1 (log(std::abs(a)) + log(std::abs(c)), dd);
     // complex_t _g1= log(a * std::conj(c));
//      INFO("g1 = "<<G1<<" dd = "<<dd<<" a = "<<a<<" c = "<<c);
      //complex_t _g1 = std::log(A.getVal(mc[i]) * std::conj(C.getVal(mc[i])));
      g1.push_back(_g1);
      G1 += _g1;
  }
  clockG1.stop();

  INFO("Got "<<g1.size()<<" g1, sum(g1) = G1 = "<<G1<<" took "<<clockG1);

  INFO("sumDD = "<<sumDD);

  ProfileClock clockG2;
  clockG2.start();
  for (size_t m=0;m<pc.getOrder()+1;m++){
    for (size_t n=0;n<pc.getOrder() + 1 - m; n++){
        std::pair<size_t, size_t> mn({m, 2*n+1});
        std::vector<double> g2_mn;
        double G2_mn=0;
        for (size_t i=0;i<mc.size();i++){
            std::vector<double> XY = pc.getXY(mc[i]);
            double _g2_mn = pc.Poly2D(XY[0], XY[1], m, 2*n+1);
            g2_mn.push_back(_g2_mn);
            G2_mn += _g2_mn;
        }

        g2.insert(std::pair<std::pair<size_t, size_t>, std::vector<double> >({mn, g2_mn})  );
        G2.insert(std::pair<std::pair<size_t, size_t>, double>({mn, G2_mn})  );
    }
  } 
  clockG2.stop();

  INFO("Got "<<g2.size()<<" took "<<clockG2);

  ProfileClock clockDelta1;
  clockDelta1.start();
  for (size_t i=0;i<mc.size();i++){
      Delta1.push_back( g1[i] - G1 );
  }
  clockDelta1.stop();
  INFO("Got "<<Delta1.size()<<" took "<<clockDelta1);

  ProfileClock clockDelta2 ;
  clockDelta2.start();
  for (size_t m=0;m<pc.getOrder()+1;m++){
      for(size_t n=0;n<pc.getOrder()+1 -m;n++){
          std::pair<size_t, size_t> mn({m,2*n+1});
          std::vector<double> Delta2_mn;
          for(size_t i=0;i<mc.size();i++){
              Delta2_mn.push_back(g2[mn][i] - G2[mn]);
          }
          Delta2.insert(std::pair<std::pair<size_t, size_t>, std::vector<double> >({mn, Delta2_mn}));
      }
  }
  clockDelta2.stop();

  INFO("Got "<<Delta2.size()<<" took "<<clockDelta2);
 
  std::vector<complex_t> Delta;
  ProfileClock clockDelta;
  clockDelta.start();
  for (size_t i=0;i<Delta1.size();i++){
    complex_t _Delta_i = Delta1[i];  
    /*
    for (auto p : Delta2){
        std::pair<size_t, size_t> ij = p.first;
        size_t m = ij.first;
        size_t n = ij.second;
        _Delta_i += complex_t(0, MPS["PhaseCorrection::C" + std::to_string(m) + "_" + std::to_string(n)]->mean() * p.second[i]);
        
    }
    */
    Delta.push_back(_Delta_i);
  }
  clockDelta.stop();
  INFO("Got "<<Delta.size()<<" Delta2 took "<<clockDelta);
  ProfileClock clockDeltaMN;
  clockDeltaMN.start();

  complex_t expDelta1; 
    for (auto p : Delta2){
        std::pair<size_t, size_t> ij = p.first;
        size_t m = ij.first;
        size_t n = ij.second;
        for (size_t i=0;i<Delta.size();i++){     
//            INFO(i<<" "<<Delta[i]<<" "<<MPS["PhaseCorrection::C" + std::to_string(m) + "_" + std::to_string(n)]->mean() * p.second[i]);
            double DMN = MPS["PhaseCorrection::C" + std::to_string(m) + "_" + std::to_string(n)]->mean() * p.second[i];
            expDelta1 += exp(Delta[i] + complex_t(0, DMN));
        }
    }

  clockDeltaMN.stop();
  INFO("DeltaMN took "<<clockDeltaMN<<" expDelta = "<<expDelta1);

  ProfileClock clockG0;
  clockG0.start();
  complex_t G0;// = G1;
  for (auto p : G2){
        std::pair<size_t, size_t> ij = p.first;
        size_t m = ij.first;
        size_t n = ij.second;
        double DMN = MPS["PhaseCorrection::C" + std::to_string(m) + "_" + std::to_string(n)]->mean() * p.second;
        G0 += (G1 + complex_t(0, DMN));
  }
  clockG0.stop();

  INFO("G0 Took "<<clockG0<<" G0 = "<<G0<<" expG0 = "<<exp(G0));
 
  complex_t inter =  exp(G0) * expDelta1 ;
  double myNorm = -2 * std::real(exp(G0) * expDelta1) + nA + nC;
 INFO("My Norm (from calculating logs and exps) = "<<myNorm<<" inter = "<<inter);
 

  return 0;
  auto G = [&MPS, &G1, &G2](){
      complex_t _G = G1;
      for (auto p : G2){
          std::pair<size_t, size_t> ij = p.first;
          size_t m = ij.first;
          size_t n = ij.second;
          _G += complex_t(0,MPS["PhaseCorrection::C" + std::to_string(m) + "_" + std::to_string(n)]->mean() * p.second);
      }
      return _G;
  };

  auto expDelta = [&MPS, &Delta1, &Delta2](){
      //std::vector<complex_t> _Delta;
      complex_t _expDelta = 0;
      for (size_t i=0;i<Delta1.size();i++){
        complex_t _Delta_i = Delta1[i];  
        for (auto p : Delta2){
            std::pair<size_t, size_t> ij = p.first;
            size_t i = ij.first;
            size_t j = ij.second;
            _Delta_i += complex_t(0, MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(j)]->mean() * p.second[i]);
            
        }
        _expDelta += std::exp(_Delta_i);
      }
      return _expDelta;
  };

  INFO("Getting G");
  ProfileClock clockG;
  clockG.start();
  complex_t G01 = G();
  clockG.stop();


  INFO("Getting Delta");
  ProfileClock clockexpDelta;
  clockDelta.start();
  complex_t expDelta0 = expDelta();
  clockexpDelta.stop();

  double myNorm2 = -2 * std::real(std::exp(G0) * expDelta0) + nA + nC;
  INFO("My Norm (from calculating logs and exps) = "<<myNorm);
  return 0;










}