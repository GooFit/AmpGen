#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/gCoherentSum.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
//#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/SumPDF.h"
//#include "AmpGen/CombCorrLL.h"
//#include "AmpGen/CombGamCorrLL.h"
//#include "AmpGen/CombLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>

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


   auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
 
  

  auto BTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();
 

  MinuitParameterSet MPS;
  MPS.loadFromStream();

  EventType eventType = EventType(pNames);

  for (auto& BTag : BTags){

INFO("B DecayType = "<<BTag);
 
 
    

    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag,' ')[1];
    int B_Conj = std::stoi(split(BTag,' ')[2]);
    int gammaSign = std::stoi(split(BTag,' ')[3]);
   bool useXY = std::stoi(split(BTag,' ')[4]);
 
    
    INFO("GammaSign = "<<gammaSign);
    if (B_Conj == 1){
      eventType = eventType.conj(true);
    }


  

  auto sig = gCoherentSum(eventType, MPS);

  std::string DataFile = NamedParameter<std::string>("DataSample", "");
  std::string IntFile = NamedParameter<std::string>("IntegrationSample", "");

  std::stringstream DataSS;
  DataSS<<DataFile<<":"<<B_Name;
  std::string DataLoc = DataSS.str();


  std::stringstream IntSS;
  IntSS<<IntFile<<":"<<B_Name;
  std::string IntLoc = IntSS.str();



  EventList Data = EventList(DataLoc, eventType);
  EventList Int = EventList(IntLoc, eventType);
  
  sig.setEvents(Data);
  sig.setMC(Int);

  auto ll = make_likelihood(Data, sig);
  Minimiser mini = Minimiser(ll, &MPS);
  mini.gradientTest();
  mini.doFit();

  std::ofstream fitOut;
  fitOut.open("fitB.csv");
  fitOut<<"x+ "<<MPS["gCoherentSum::x"]->mean()<<" "<<MPS["gCoherentSum::x"]->err()<<"\n"<<"y+ "<<MPS["gCoherentSum::y"]->mean()<<" "<<MPS["gCoherentSum::y"]->err()<<"\n";
  fitOut.close();

  }


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


