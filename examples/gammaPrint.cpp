#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
//#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/SumPDF.h"
//#include "AmpGen/CombCorrLL.h"
//#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/CombLL.h"
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

  std::vector<EventList> SigData;
  std::vector<EventList> SigInt;
  std::vector<EventType> SigType;
  std::vector<std::string> sumFactors;
  std::vector<int> gammaSigns;
  std::vector<int> useXYs;
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

    auto sig = pCoherentSum(eventType, MPS ,B_Pref, gammaSign, useXY);



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


    SigData.push_back(Data);
    SigInt.push_back(Int);
    SigType.push_back(eventType);
    sumFactors.push_back(B_Pref);
    gammaSigns.push_back(gammaSign);
    useXYs.push_back(useXY);

    sig.prepare();
  std::vector<double> ci = {};
  std::vector<double> si = {};
  std::vector<double> bins = {};
  for (int j=1;j<9;j++){
    double _ci =0;
    double _si = 0;
    double _N = 0;
    auto eLow = 2 * M_PI * j/8. - 5 * M_PI/4.;
    auto eHigh = 2 * M_PI * (j+1)/8. - 5 * M_PI/4.;
    std::ofstream outBinned;
    std::ofstream outBinnedM;

    std::stringstream fileNamecisi;
    std::stringstream fileNameBinned;
    std::stringstream fileNameBinnedM;
    fileNameBinned<<"Bin"<<B_Name<<j<<".csv";
    fileNameBinnedM<<"Bin"<<B_Name<<-j<<".csv";

    outBinned.open(fileNameBinned.str());
    outBinnedM.open(fileNameBinnedM.str());

    for (int i=0;i<Data.size();i++){
      auto event = Data[i];
      auto x = event.s(0,1);
      auto y = event.s(0,2);
      auto z = event.s(1,2);
      auto mp = sqrt(event.s(1,1))/2.;
      auto mm = sqrt(event.s(2,2))/2.;
      auto mK = sqrt(event.s(0,0))/2.;
      auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
      Expression xmin = pow(mp + mK, 2);
      Expression xmax = pow(mD - mm, 2);
      Expression x0 = (xmax + xmin)/2;
      Expression ymin = pow(mp + mK, 2);
      Expression ymax = pow(mD - mp, 2);
      Expression y0 = (ymax + ymin)/2;
      Expression X = (2 * x - xmax - xmin)/(xmax - xmin);
      Expression Y = (2 * y - ymax - ymin)/(ymax - ymin);
      auto A = sig.A(event);
      auto C = sig.C(event);
      auto psi = sig.getVal(event);
      auto corr = sig.correction(event);
      auto ACstCorr = A * std::conj(C) * fcn::exp(Constant(0, 1) * corr)();

      auto dd = std::arg(ACstCorr);

      if (dd > eLow && dd < eHigh){
        if (x>y){
          outBinned << x <<"\t"
          << y << "\t"
          << dd << "\t"
          << psi.real()<<"\t"
          << psi.imag()<<"\t"
          << A.real()<<"\t"
          << A.imag()<<"\t"
          << C.real()<<"\t"
          << C.imag()<<"\t"
          <<"\n";          

          _ci += fcn::cos(dd)().real();
          _si += fcn::sin(dd)().real();
          _N +=1;

        }
        else{
          outBinnedM << x <<"\t"
          << y << "\t"
          << dd << "\t"
          << psi.real()<<"\t"
          << psi.imag()<<"\t"
          << A.real()<<"\t"
          << A.imag()<<"\t"
          << C.real()<<"\t"
          << C.imag()<<"\t"
          <<"\n";          

          //_ci += fcn::cos(dd)().real();
          //_si += fcn::sin(dd)().real();
          //_N +=1;


        }

          INFO(x <<" "
          << y << " "
          << dd << " "
          << psi.real()<<" "
          << psi.imag()<<" "
          << A.real()<<" "
          << A.imag()<<" "
          << C.real()<<" "
          << C.imag()<<" "
          <<j);
      }
    }
    ci.push_back(_ci/_N);
    si.push_back(_si/_N);
    bins.push_back(j);
    outBinned.close();
    outBinnedM.close();
  }
    std::ofstream outcisi;
    outcisi.open("Bcisi.csv");
    for (int i=0;i<bins.size();i++){
      outcisi<<bins[i]<<"\t"<<ci[i]<<"\t"<<si[i]<<"\n";
    }
    outcisi.close();
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


