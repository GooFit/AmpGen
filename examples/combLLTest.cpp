#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/CombCorrLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>


//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps){
  std::map<std::string, std::vector<double> > out;
  for (auto& param : mps){
    std::vector<double> vect = {param->mean(), param->err()};
    out[param->name()] = vect;
  }
  return out;
}

MinuitParameterSet * copyMPS(MinuitParameterSet& mps){
  MinuitParameterSet * newMPS = new MinuitParameterSet; 
  for (auto& param : mps){
    newMPS->add(param);
  }
  return newMPS;
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

void add_CP_conjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  for( auto& param : mps ){
    const std::string name = param->name();
    size_t pos=0;
    std::string new_name = name; 
    int sgn=1;
    Flag flag = param->flag();
    //Flag flag = Flag::Fix;
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


int main( int argc, char* argv[] )
{
      OptionsParser::setArgs( argc, argv );

    size_t hwt = std::thread::hardware_concurrency();
    size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
    
    
    size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
    
    
    
    std::string output  = NamedParameter<std::string>("Output" , "ToyMC.root", "File containing output events"); 
    auto pNames = NamedParameter<std::string>("EventType" , ""    
        , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
    auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
    bool m_debug        = NamedParameter<bool>("Debug", false, "Debug QcFitter output");
    
    int nBins = NamedParameter<int>("nBins", 100, "number of bins for projection");
    int nFits = NamedParameter<int>("nFits", 4, "number of repeats of mini.doFits() for debug purposes!");
    bool doProjections = NamedParameter<bool>("doProjections", true);
    bool doPCorrSum = NamedParameter<bool>("doPCorrSum", false);
    std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");

    bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
    bool doCombFit      = NamedParameter<bool>("doCombFit", false, "Do combined fit");
    bool doTagFit      = NamedParameter<bool>("doTagFit", false, "Do fit for each tag");
    bool QcGen2 = NamedParameter<bool>("QcGen2", false, "internal boolean - for new QcGenerator");
    bool doFit = NamedParameter<bool>("doFit", true, "Do the fit");
    if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
    if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);
    if (intFile == ""){

    }

      TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  if (m_debug) INFO("LogFile: " << logFile << "; Plots: " << plotFile );
   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  if (m_debug) INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
    #endif

  std::vector<std::string> varNames = {"E", "PX", "PY", "PZ"};
  //auto yc = DTYieldCalculator(crossSection);
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(MPS);
  }
  std::map<std::string, std::vector<double> > inits = getParams(MPS);


    INFO("Doing loop of Fits");

    std::vector<EventList> SigData;
    std::vector<EventList> TagData;
    std::vector<EventList> SigInt;
    std::vector<EventList> TagInt;
    std::vector<EventType> SigType;
    std::vector<EventType> TagType;
    std::vector<std::string> sumFactors;
    double LL_comb_sum = 0;
    std::vector<double> LL_tags;
    for (int i=0; i < tags.size(); i++){
    MinuitParameterSet * MPS_tag = new MinuitParameterSet();
    MPS_tag->loadFromStream();
    if (makeCPConj){
        INFO("Making CP conjugate states");
        add_CP_conjugate(*MPS_tag);
    }
        std::stringstream tag_log;
        auto tagName = split(tags[i],' ')[0];
        tag_log<<tagName<<"_fit.log";
        std::stringstream tag_fit;
        tag_fit<<tagName<<"_plots.root";
        auto tag_plotName = tag_fit.str();
        auto tag_logName = tag_log.str();

        auto sigevents_tag = getEvents("signal", pNames, tags[i], dataFile, intFile);
        auto sigMCevents_tag = getEvents("sigMC", pNames, tags[i], dataFile, intFile);
        auto tagevents_tag = getEvents("tag", pNames, tags[i], dataFile, intFile);
        auto tagMCevents_tag = getEvents("tagMC", pNames, tags[i], dataFile, intFile);

        SigData.push_back(sigevents_tag);
        TagData.push_back(tagevents_tag);
        SigInt.push_back(sigMCevents_tag);
        TagInt.push_back(tagMCevents_tag);
        SigType.push_back(sigevents_tag.eventType());
        TagType.push_back(tagevents_tag.eventType());
        sumFactors.push_back("Psi3770");

      auto cs_tag = pCorrelatedSum(sigevents_tag.eventType(), tagevents_tag.eventType(), *MPS_tag);
      cs_tag.setEvents(sigevents_tag, tagevents_tag);
      cs_tag.setMC(sigMCevents_tag, tagMCevents_tag);
      cs_tag.prepare();
      //auto LL_tag2 = make_likelihood( events_tag["signal"], events_tag["tag"], false, cs_tag);
      auto LL_tag2 = make_likelihood( sigevents_tag, tagevents_tag, cs_tag);
      auto mini_tag = Minimiser(LL_tag2, MPS_tag);
      mini_tag.prepare();
      auto LL_tag = LL_tag2.getVal();
      LL_comb_sum += LL_tag;
      LL_tags.push_back(LL_tag);
      INFO("LL for "<<tagName<<" = "<<LL_tag);
    }

    CombCorrLL combLL = CombCorrLL(SigData, TagData, SigInt, TagInt, SigType, TagType, MPS, sumFactors);
    INFO("Making Combined Minimiser object");
    Minimiser combMini = Minimiser(combLL, &MPS);
    combMini.prepare();
    INFO("Minimising now");
    double LL_comb = combLL.getVal();
    for (auto LL : LL_tags){
        INFO("LL = "<<LL);
    }
    INFO("LL_comb = "<< LL_comb);
    INFO("Sum of tag LL = "<<LL_comb_sum);
    INFO("Difference = "<<LL_comb - LL_comb_sum);
    return 0;
}