#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
//#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/SumPDF.h"
//#include "AmpGen/CombCorrLL.h"
//#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/CombLL.h"
#include "AmpGen/GamLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>
#include "AmpGen/AddCPConjugate.h"
//#include <boost/algorithm/string.hpp>
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


  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */


  bool m_debug = NamedParameter<bool>("debug", false);
  int nThreads = NamedParameter<int>("nThreads", 12);

   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  if (m_debug) INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
  #endif

  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;
   auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
 
  

  auto BTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();
 

  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj)       AddCPConjugate(MPS);

  EventType eventType = EventType(pNames);

  std::vector<EventList> SigData;
  std::vector<EventList> SigInt;
  std::vector<EventType> SigType;
  std::vector<std::string> sumFactors;
  std::vector<int> gammaSigns;
  std::vector<int> useXYs;
  std::vector<int> B_Conjs;
  bool fitEach = NamedParameter<bool>("FitEach", true);
  SimFit simfit;
  std::vector<SumPDF<EventList, pCoherentSum&>> pdfs;

int NInt = NamedParameter<int>("NInt", 1e7);

    
    INFO("Using "<<NInt<<" integration events");

//  EventList mc =  Generator<>(eventType, &rndm).generate(NInt);

  for (auto& BTag : BTags){

    INFO("B DecayType = "<<BTag);
    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag,' ')[1];
    int B_Conj = std::stoi(split(BTag,' ')[2]);
    int gammaSign = std::stoi(split(BTag,' ')[3]);
    bool useXY = std::stoi(split(BTag,' ')[4]);
 
    
    INFO("GammaSign = "<<gammaSign);
    if (B_Conj == 1){
      //eventType = eventType.conj();
    }



    std::string DataFile = NamedParameter<std::string>("DataSample", "");
    std::string IntFile = NamedParameter<std::string>("IntegrationSample", "");

    std::stringstream DataSS;
    DataSS<<DataFile<<":"<<B_Name;
    std::string DataLoc = DataSS.str();


    std::stringstream IntSS;
    IntSS<<IntFile<<":"<<B_Name;
    std::string IntLoc = IntSS.str();

    EventList Data = EventList(DataLoc, eventType);
//    EventList Int = EventList(IntLoc, eventType);

    
    SigData.emplace_back(Data);
    //SigInt.emplace_back(mc);
    SigType.emplace_back(eventType);
    sumFactors.emplace_back(B_Pref);
    gammaSigns.emplace_back(gammaSign);
    useXYs.emplace_back(useXY);
    B_Conjs.emplace_back(B_Conj);


  }

 //auto LLC = CombLL(SigData, mc, SigType, MPS, sumFactors, gammaSigns, useXYs, B_Conjs);
//
/*
  pdfs.reserve(SigData.size());
  std::vector<pCoherentSum> fcs(SigData.size());
  for (size_t i=0;i<SigData.size(); i++){
   fcs[i] = pCoherentSum(eventType, MPS, sumFactors[i], gammaSigns[i], useXYs[i], B_Conjs[i]);
   pdfs.emplace_back( make_pdf(fcs[i])); 
   pdfs[i].setEvents(SigData[i]);
//   auto& mc  = SigInt[i];
   for_each(pdfs[i].pdfs(), [&mc](auto& pdf){pdf.setMC(mc);});

   simfit.add(pdfs[i]);


    }
*/   
/*
  Minimiser mini(LLC, &MPS);
mini.prepare();
INFO("Mini = "<<mini.FCN());
  mini.gradientTest();
  mini.doFit();
  FitResult * fr = new FitResult(mini);
  fr->writeToFile(NamedParameter<std::string>("Logfile", "BFit.log"));
  auto covMatrix = mini.covMatrix();
  std::string covMatrixFile = NamedParameter<std::string>("CovOutput", "BCov.root");
  TFile * fCov = new TFile(covMatrixFile.c_str(), "recreate");
  fCov->cd();
  covMatrix.Write("CovMatrix");
  delete fCov;
  */

  GamLL LL (SigData, eventType, MPS, gammaSigns, useXYs, B_Conjs);
  auto sf0 = LL.sumFactor(gammaSigns[0], useXYs[0]);
  auto sf1 = LL.sumFactor(gammaSigns[1], useXYs[1]);

  INFO("sf0 = "<<sf0);
  INFO("sf1 = "<<sf1);

  auto A = LL.get_A();
  auto AMC = LL.get_AMC();

  real_t n0 = LL.norm(0);
  real_t n1 = LL.norm(1);
  INFO("norm0 = "<<n0);
  INFO("norm1 = "<<n1);
  real_t ll = LL();
  INFO("ll = "<<ll);
  return 0;
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


