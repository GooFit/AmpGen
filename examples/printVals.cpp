#include <Rtypes.h>
#include <TH1.h>
#include <dlfcn.h>
#include <memory>
#include <string>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "TGraph2D.h"



#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#include "AmpGen/DynamicFCN.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/RecursivePhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/QcGenerator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/enum.h"

#include "AmpGen/Psi3770.h"


#include "AmpGen/MetaUtils.h"
#include "AmpGen/corrEventList.h"

#include "AmpGen/Kinematics.h"


using namespace AmpGen;
namespace AmpGen { make_enum(generatorType, CoherentSum, PolarisedSum, FixedLib, RGenerator) }

void add_CP_conjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  for( auto& param : mps ){
    const std::string name = param->name();
    size_t pos=0;
    std::string new_name = name; 
    int sgn=1;
    Flag flag = param->flag();
    if( name.find("::") != std::string::npos ){
      pos = name.find("::");
      auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
      if( props != 0 ) new_name = props->anti().name() + name.substr(pos); 
    }
    else { 
      auto tokens=split(name,'_');
      std::string reOrIm = *tokens.rbegin();
      std::string name   = tokens[0];

      if ( reOrIm == "Re" || reOrIm == "Im" ){
        auto p = Particle( name ).conj();
        sgn = reOrIm == "Re" ? p.CP() : 1; 
        new_name = p.uniqueString() +"_"+reOrIm;
      }
      else if( tokens.size() == 2 ) {
        auto props = AmpGen::ParticlePropertiesList::get( name );
        if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
      }
    }
    if( mps.find( new_name ) == nullptr ){
      tmp.push_back( new MinuitParameter(new_name, flag, sgn * param->mean(), param->err(), 0, 0));
    }
  }
  for( auto& p : tmp ) mps.add( p );
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




std::map<std::string, EventType> makeEventTypes(std::vector<std::string> sigName, std::string tagName){
  EventType signalType( sigName );
  INFO("Getting Tag name split "<<tagName);
  auto tokens       = split(tagName, ' ');
  INFO("Building Tag particle"<<tokens[1]);
  auto tagParticle  = Particle(tokens[1], {}, false);
  INFO("Build EventType");
  EventType    tagType = tagParticle.eventType(); 
  auto sigBranches = makeBranches(signalType, "");
  auto tagBranches = makeBranches(tagType, "Tag_");

  std::map<std::string, EventType> types = {};

  types["signal"] = signalType;
  types["tag"] = tagType;
  return types;

}


  
EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile, std::string intFile){
   
    /*
    EventType signalType( sigName );
    INFO("Getting Tag name split "<<tagName);
    auto tokens       = split(tagName, ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto tagParticle  = Particle(tokens[1], {}, false);
    INFO("Build EventType");
    EventType    tagType = tagParticle.eventType(); 
    */

    INFO("Getting Tag name split "<<tagName);
    auto tokens       = split(tagName, ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto types = makeEventTypes(sigName, tagName);
    INFO("Build EventType");
    auto signalType = types["signal"];
    auto tagType = types["tag"];
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");

    EventList sigEvents;
    EventList tagEvents;
    EventList sigMCEvents;
    EventList tagMCEvents;
    EventList null;


    INFO("Generated Events used QcGen2");
    std::stringstream signame;
    signame<<":Signal_";
    signame<<tokens[0];
    std::stringstream tagname;
    tagname<<":Tag_";
    tagname<<tokens[0];
    sigEvents = EventList(dataFile + signame.str() , signalType);
    tagEvents = EventList(dataFile + tagname.str() , tagType);
    sigMCEvents = EventList(intFile + signame.str() , signalType);
    tagMCEvents = EventList(intFile + tagname.str() , tagType);
    
    if (type=="signal"){
     return sigEvents;
    }
    else if (type=="tag"){ 
        return tagEvents;
    }
    else if (type=="sigMC"){
        return sigMCEvents;
    }
    else if (type=="tagMC") {
        return tagMCEvents;
    }
    else {
        return null;
    }
    


}



int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 100, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );


   bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");

  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    INFO("Using: " << nCores  << " / " << concurentThreadsSupported  << " threads" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif 
  EventType sigType( pNames );

  bool outputVals = NamedParameter<bool>("outputVals", false, "Whether to print out values for amplitude in a csv file");
  std::string ampFile = NamedParameter<std::string>("ampFile", "corr.csv", "Output file for correlated amplitude");
  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );
  bool doBoost  = NamedParameter<bool>("doBoost", true, "Boost to psi(3770) frame");
  bool doQC  = NamedParameter<bool>("doQC", true, "Boost to psi(3770) frame");
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();


  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile = NamedParameter<std::string>("IntegrationSample", ""          , "Name of file containing flat sample." );




  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

if (doQC){
for (int i=0; i < tags.size(); i++){
   MinuitParameterSet *  MPS = new MinuitParameterSet(); 
  MPS->loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(*MPS);
  }
    auto sigevents_tag = getEvents("signal", pNames, tags[i], dataFile, intFile);
    auto sigMCevents_tag = getEvents("sigMC", pNames, tags[i], dataFile, intFile);
    auto tagevents_tag = getEvents("tag", pNames, tags[i], dataFile, intFile);
    auto tagMCevents_tag = getEvents("tagMC", pNames, tags[i], dataFile, intFile);

    auto sigEventType = sigevents_tag.eventType();
    auto tagEventType = tagevents_tag.eventType();
    CoherentSum sig(sigEventType, *MPS);
    pCorrelatedSum cs_tag (sigevents_tag.eventType(), tagevents_tag.eventType(), *MPS);

     
    

    cs_tag.setEvents(sigevents_tag, tagevents_tag);
    cs_tag.setMC(sigMCevents_tag, tagMCevents_tag);
    cs_tag.prepare();
 std::ofstream out;



 
    
    std::stringstream tag_log;
    auto tagName = split(tags[i],' ')[0];
    tag_log<<tagName<<"_vals.csv"; 
    auto tag_ampFile = tag_log.str();
    out.open(tag_ampFile.c_str());
    auto norm = cs_tag.norm();
    auto aA = cs_tag.getA();
    auto aB = cs_tag.getB();
    auto aC = cs_tag.getC();
    auto aD = cs_tag.getD();

    double slowNorm = 0;
    double testA = 0;
    double testC = 0;
    for (int i=0;i<sigMCevents_tag.size();i++){
      slowNorm += std::norm(cs_tag.getVal(sigMCevents_tag[i], tagMCevents_tag[i]))/norm;
      testA += aA(sigMCevents_tag[i])/sigMCevents_tag.size();
      testC += aC(sigMCevents_tag[i])/sigMCevents_tag.size();
    }
    double testData = 0;
    double testDataA = 0;
    double testDataC = 0;

    for (int i=0; i < sigevents_tag.size(); i++){
    auto evt_sig = sigevents_tag[i];
    auto evt_tag = tagevents_tag[i];
    testData += std::norm(cs_tag.getVal(evt_sig, evt_tag));
    testDataA += std::norm(aA.getVal(evt_sig));
    testDataC += std::norm(aC.getVal(evt_sig));
    }
INFO("testData = "<<testData);
INFO("testDataA = "<<testDataA);
INFO("testDataC = "<<testDataC);
    INFO("Norm = "<<norm);
    INFO("ANorm = "<<aA.norm());
    INFO("CNorm = "<<aC.norm());
    INFO("testA = "<<testA);
    INFO("testC = "<<testC);
    INFO("slowNorm = "<<slowNorm/sigMCevents_tag.size());
    INFO("normData? = "<<testData*norm/sigevents_tag.size());
    double normDat = (double)sigMCevents_tag.size()/(double)sigevents_tag.size();

    INFO("normData = "<<normDat);
        auto rnorm = pow(norm, 0.5);
    for (int i=0; i < sigevents_tag.size(); i++){

    auto evt_sig = sigevents_tag[i];
    auto evt_tag = tagevents_tag[i];
    testData += cs_tag.prob(evt_sig, evt_tag);
    auto s01 = evt_sig.s(0,1);
    auto s02 = evt_sig.s(0,2);
    auto s12 = evt_sig.s(1,2);
    auto vals = cs_tag.getVals(evt_sig, evt_tag);
    auto A = aA.getVal(evt_sig);
    auto B = aB.getVal(evt_tag);
    auto C = aC.getVal(evt_sig);
    auto D = aD.getVal(evt_tag);

  auto A2 = aA(evt_sig)/sigevents_tag.size();
    auto B2 = aB(evt_tag)/sigevents_tag.size();
    auto C2 = aC(evt_sig)/sigevents_tag.size();
    auto D2 = aD(evt_tag)/sigevents_tag.size();


    auto ABCD = cs_tag.getVal(evt_sig, evt_tag);
    auto corr = vals[5];
    auto ACst = A * std::conj(C)/norm;
    auto dCorr = cs_tag.errcorrection(evt_sig);

    auto prob = cs_tag.prob(evt_sig, evt_tag)/sigevents_tag.size();
    auto cosDD = -((fcn::pow(abs(ABCD), 2)() - fcn::pow(abs(A), 2)() * fcn::pow(abs(B), 2)() - fcn::pow(abs(C), 2)() * fcn::pow(abs(D), 2)())/(2 * abs(A) * abs(B) * abs(C) * abs(D) )).real();

    auto dd = std::arg(ACst);
    auto ACstCorr = A * std::conj(C) * fcn::exp(Constant(0, 1) * corr)();
    auto ddCorr = std::arg(ACstCorr);
    
    out << s01 << "\t" //0
        << s02 << "\t" //1
        << s12 << "\t" //2
        << ddCorr  << "\t" //3
        << cosDD << "\t" //4
        <<corr.real()<<"\t" //5
        <<dCorr.real()<<"\t" //6
        << dd<<"\t" //7
        <<std::norm(ABCD)<<"\t" //8
        <<A.real()<<"\t" //9
        <<A.imag()<<"\t" //10 
        <<std::abs(A)<<"\t" //11
        <<std::arg(A)<<"\t" //12
        <<C.real()<<"\t" //13
        <<C.imag()<<"\t" //14
        <<std::abs(C)<<"\t" //15
        <<std::arg(C)<<"\t" //16

        <<ABCD.real()<<"\t"
        <<ABCD.imag()<<"\t"

        <<prob<<"\t"
        <<A2<<"\t"
        <<B2<<"\t"
        <<C2<<"\t"
        <<D2<<"\t"

        <<"\n";

      //Binning
}
  out.close();

  


for (int j=1;j<9;j++){
  auto eLow = 2 * M_PI * j/8. - 5 * M_PI/4.;
  auto eHigh = 2 * M_PI * (j+1)/8. - 5 * M_PI/4.;
  std::ofstream outBinned;
  std::ofstream outBinnedM;
  std::ofstream outcisi;
  std::stringstream fileNamecisi;
  std::stringstream fileNameBinned;
  std::stringstream fileNameBinnedM;
  fileNameBinned<<"Bin"<<tagName<<j<<".csv";
  fileNameBinnedM<<"Bin"<<tagName<<-j<<".csv";
  fileNamecisi<<"cisi"<<tagName<<j<<".csv";
     auto aA = cs_tag.getA();
    auto aB = cs_tag.getB();
    auto aC = cs_tag.getC();
    auto aD = cs_tag.getD();
 
  outBinned.open(fileNameBinned.str());
  outBinnedM.open(fileNameBinnedM.str());
  outcisi.open(fileNamecisi.str());
  auto ci = 0.0;
  auto si =0.0 ;
  auto Np = 0.0;
  auto Nm = 0.0;
  for (int i=0; i < sigevents_tag.size(); i++){

    auto evt_sig = sigevents_tag[i];
    auto evt_tag = tagevents_tag[i];
    auto s01 = evt_sig.s(0,1);
    auto s02 = evt_sig.s(0,2);
    auto s12 = evt_sig.s(1,2);
    auto vals = cs_tag.getVals(evt_sig, evt_tag);

    auto A = aA.getVal(evt_sig);
    auto B = aB.getVal(evt_tag);
    auto C = aC.getVal(evt_sig);
    auto D = aD.getVal(evt_tag);


    auto ABCD = cs_tag.getVal(evt_sig, evt_tag);
    auto corr = vals[5];
    auto prob = cs_tag.prob(evt_sig, evt_tag);
//    INFO("prob = "<<prob<<" Aprob = "<<aA(evt_sig)<<" Cprob = "<<aC(evt_sig));
    auto ACst = A * std::conj(C)/norm;
    auto dCorr = cs_tag.errcorrection(evt_sig);
  auto A2 = aA(evt_sig);
    auto B2 = aB(evt_tag);
    auto C2 = aC(evt_sig);
    auto D2 = aD(evt_tag);



    auto cosDD = -((fcn::pow(abs(ABCD), 2)() - fcn::pow(abs(A), 2)() * fcn::pow(abs(B), 2)() - fcn::pow(abs(C), 2)() * fcn::pow(abs(D), 2)())/(2 * abs(A) * abs(B) * abs(C) * abs(D) )).real();

    auto dd = std::arg(ACst);
    auto ACstCorr = A * std::conj(C) * fcn::exp(Constant(0, 1) * corr)();
    auto ddCorr = std::arg(ACstCorr);
 
if (dd > eLow && dd < eHigh){

if (s01>s02){
    outBinned << s01 << "\t" //0
        << s02 << "\t" //1
        << s12 << "\t" //2
        << ddCorr  << "\t" //3
        << cosDD << "\t" //4
        <<corr.real()<<"\t" //5
        <<dCorr.real()<<"\t" //6
        << dd<<"\t" //7
        <<std::norm(ABCD)<<"\t" //8
        <<A.real()<<"\t" //9
        <<A.imag()<<"\t" //10 
        <<std::abs(A)<<"\t" //11
        <<std::arg(A)<<"\t"
        <<C.real()<<"\t"
        <<C.imag()<<"\t"
        <<std::abs(C)<<"\t"
        <<std::arg(C)<<"\t"
        <<ABCD.real()<<"\t"
        <<ABCD.imag()<<"\t"
        <<prob<<"\t"
        <<A2<<"\t"
        <<B2<<"\t"
        <<C2<<"\t"
        <<D2<<"\t"



        <<"\n";

        ci += fcn::cos(dd)().real();

       
        si += fcn::sin(dd)().real();
        
        Np+=1;
}
else{
    outBinnedM << s01 << "\t" //0
        << s02 << "\t" //1
        << s12 << "\t" //2
        << ddCorr  << "\t" //3
        << cosDD << "\t" //4
        <<corr.real()<<"\t" //5
        <<dCorr.real()<<"\t" //6
        << dd<<"\t" //7
        <<std::norm(ABCD)<<"\t" //8
        <<A.real()<<"\t" //9
        <<A.imag()<<"\t" //10 
        <<std::abs(A)<<"\t" //11
        <<std::arg(A)<<"\t"
        <<C.real()<<"\t"
        <<C.imag()<<"\t"
        <<std::abs(C)<<"\t"
        <<std::arg(C)<<"\t"
        <<ABCD.real()<<"\t"
        <<ABCD.imag()<<"\t"

        <<prob<<"\t"
        <<A2<<"\t"
        <<B2<<"\t"
        <<C2<<"\t"
        <<D2<<"\t"



        <<"\n";
        Nm +=1;
        ci += fcn::cos(dd)().real();

       
        si += -fcn::sin(dd)().real();
 

}



}



}
INFO("N+ = "<<Np);
INFO("N- = "<<Nm);
INFO("N = "<<Np + Nm);

outcisi<<ci/(Np+Nm)<<"\t"<<si/(Np+Nm)<<"\n";
outcisi.close();
outBinned.close();
outBinnedM.close();



}
}
}
else{
   MinuitParameterSet *  MPS = new MinuitParameterSet(); 
  MPS->loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(*MPS);
  }
   EventList events(dataFile, eventType);

  EventList eventsMC =  Generator<>(eventType, &rndm).generate(2e6);
  CoherentSum sig(eventType, *MPS);
  sig.setEvents(events);
  sig.setMC( eventsMC );
  sig.prepare();
     std::ofstream out;
    out.open(ampFile.c_str());
 
  auto i=0;
  for (auto& evt : events){
    auto s01 = evt.s(0,1);
    auto s02 = evt.s(0,2);
    auto s12 = evt.s(1,2);
    auto A = sig.getVal(evt);
    INFO("Event "<<i<<" s01 = "<<s01 <<" s02 = "<<s02<<" s12 = "<<s12);
    INFO("A(s01, s02) = "<<A);
    out << s01 << "\t" << s02 << "\t" << s12 << "\t" << A.real() << "\t" << A.imag()<<"\n";
    i++;

  }
    out.close();
}
return 0;
}
