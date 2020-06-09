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
    for (int i=0; i < sigevents_tag.size(); i++){

    auto evt_sig = sigevents_tag[i];
    auto evt_tag = tagevents_tag[i];
    auto s01 = evt_sig.s(0,1);
    auto s02 = evt_sig.s(0,2);
    auto s12 = evt_sig.s(1,2);
    auto vals = cs_tag.getVals(evt_sig, evt_tag);
    auto A = vals[0];
    auto B = vals[1];
    auto C = vals[2];
    auto D = vals[3];
    auto ABCD = vals[4];
    auto corr = vals[5];
    auto ACst = A * std::conj(C);
    auto dCorr = cs_tag.errcorrection(evt_sig);

    auto cosDD = -((fcn::pow(abs(ABCD), 2)() - fcn::pow(abs(A), 2)() * fcn::pow(abs(B), 2)() - fcn::pow(abs(C), 2)() * fcn::pow(abs(D), 2)())/(2 * abs(A) * abs(B) * abs(C) * abs(D) )).real();

    auto dd = std::arg(ACst);
    
    out << s01 << "\t" //0
        << s02 << "\t" //1
        << s12 << "\t" //2
        << dd  << "\t" //3
        << cosDD << "\t" //4
        <<corr.real()<<"\t" //5
        <<dCorr.real()<<"\t" //6
        << std::arg(C) - std::arg(A)<<"\t" //7
        <<std::norm(ABCD)<<"\t" //8
        <<A.real()<<"\t" //9
        <<A.imag()<<"\t" //10 
        <<std::abs(A)<<"\t" //11
        <<std::arg(A)<<"\t"
        <<C.real()<<"\t"
        <<C.imag()<<"\t"
        <<std::abs(C)<<"\t"
        <<std::arg(C)<<"\t"

        <<"\n";

}
  out.close();

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
