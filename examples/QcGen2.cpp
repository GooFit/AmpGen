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
#include "AmpGen/AddCPConjugate.h"
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

struct FixedLibPDF {
  void* lib = {nullptr};
  AmpGen::DynamicFCN<double( const double*, int )> PDF;

  void prepare(){};
  void setEvents( AmpGen::EventList& evts ){};
  double prob_unnormalised( const AmpGen::Event& evt ) const { return PDF( evt, 1 ); }
  FixedLibPDF( const std::string& lib )
  {
    void* handle = dlopen( lib.c_str(), RTLD_NOW );
    if ( handle == nullptr ) ERROR( dlerror() );
    PDF = AmpGen::DynamicFCN<double( const double*, int )>( handle, "FCN" );
  }
  size_t size() { return 0; }
  void reset( const bool& flag = false ){};
};

template <class PDF_TYPE, class PRIOR_TYPE> 
  void GenerateEvents( EventList& eventsSig
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

template <class PDF_TYPE, class PRIOR_TYPE> 
  void FilterEvents( EventList& eventsSig
                       , EventList& eventsTag
                       , EventList& inSig
                       , EventList& inTag
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
  signalGenerator.filterEventList( pdf, eventsSig, eventsTag,inSig, inTag, nEvents );
}



void add_CP_conjugate( MinuitParameterSet& mps );

EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string intFile);
std::map<std::string, EventType> makeEventTypes(std::vector<std::string> sigName, std::string tagName);
std::vector<std::string> makeBranches(EventType Type, std::string prefix);
int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 100, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
  std::string intfile = NamedParameter<std::string>("IntFile"   , "flat.root" , "Name of flat file" ); 
  std::string rwfile = NamedParameter<std::string>("RWFile"   , "RW.root" , "Name of reweighted file" ); 
  bool qcWeight	      = NamedParameter<bool>("qcWeight", false, "Do QC weighting");
  bool debug          = NamedParameter<bool>("debug", false);
  auto genType        = NamedParameter<generatorType>( "Type", generatorType::CoherentSum, optionalHelpString("Generator configuration to use:", 
    { {"CoherentSum" , "Full phase-space generator with (pseudo)scalar amplitude"}
    , {"PolarisedSum", "Full phase-space generator with particles carrying spin in the initial/final states"}
    , {"FixedLib"    , "Full phase-space generator with an amplitude from a precompiled library"}
    , {"RGenerator"  , "Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons"} } ) );
   bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");

  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
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
  std::string inputSignal = NamedParameter<std::string> ("inputSignal", "inputSignal.root", "The initial signal events");
  std::string inputTagPref = NamedParameter<std::string> ("inputTagPref", "inputTag_", "The file format for the input tag events - e.g. inputTag_ + KK = inputTag_KK.root");
  std::string sumFactor = NamedParameter<std::string>("sumFactor", "Psi3770", "KeyWord for SumFactor between A(X->f) and A(Xbar -> f) - default is Psi3770 -> A - Abar");


TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );

for( auto& tag : tags ){


    
    INFO("tag = "<<tag);
    auto tokens       = split(tag, ' ');
    INFO("tag = "<<tokens[0]);
    int nEvents_tag = nEvents;
    if (tokens.size() > 2){
      double fraction = std::stod(tokens[2]);
      nEvents_tag *= fraction;
      INFO("Generate "<<nEvents_tag);
    }
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType();
  TRandom3 rand;
  rand.SetSeed( seed + 934534 );
  //auto x = MinuitParameter("QcGen::X", Flag::Fix, 4, 0);
  //INFO("X = "<<x.mean());
  MinuitParameterSet* MPS = new MinuitParameterSet();
  MPS->loadFromStream();
  for (auto param : *MPS){
    INFO(param->name()<<" = "<<param->mean());
  }
  INFO("makeCPConj = "<<makeCPConj);
  if (makeCPConj){
    INFO("Making CP conjugate states");
    //add_CP_conjugate(*MPS);
    AddCPConjugate(*MPS);
  }

  Particle p;

  

//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );






//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );



  INFO("Generating time-dependence? " << eventType.isTimeDependent() );
  EventList acceptedSig( sigType );
  EventList acceptedTag( tagType );

  INFO("Generating events with type = " << sigType << " and "<<tagType);


 // if ( genType == generatorType::CoherentSum ) 
  INFO("Making CorrelatedSum");
    pCorrelatedSum cs( sigType, tagType, *MPS ,sumFactor);

    PhaseSpace phspSig(sigType,&rand);
    PhaseSpace phspTag(tagType,&rand);
    INFO("Generating Events now!");
    GenerateEvents( acceptedSig, acceptedTag, cs, phspSig, phspTag , nEvents_tag, blockSize, &rand );

    if (debug){
	    for (int i=0; i<acceptedSig.size(); i++){
		    INFO("Value = "<<cs.getVal(acceptedSig[i], acceptedTag[i]));
	    }
    }


      auto headPhsp = EventType({"psi(3770)0","D0","Dbar0"});
      auto psi_q = PhaseSpace(headPhsp); 
      auto beta  = [](const Event& event, const size_t&j){ return sqrt( event[4*j+0]*event[4*j+0] + event[4*j+1]*event[4*j+1] + event[4*j+2]*event[4*j+2] )/event[4*j+3] ; };
      auto p4     = [](const Event& event, const size_t&j){ return std::make_tuple(event[4*j+0], event[4*j+1], event[4*j+2]); };
      if (doBoost){
    for (size_t i=0 ; i < acceptedSig.size(); i++){
        auto psi_event = psi_q.makeEvent();
        boost( acceptedSig[i], p4(psi_event,0), beta(psi_event,0));
        boost( acceptedTag[i]   , p4(psi_event,1), beta(psi_event,1));
      }
      }

    if( acceptedSig.size() == 0 ) return -1;

    INFO( "Writing output file " );
    std::stringstream signame;
    signame<<"Signal_";
    signame<<tokens[0];
    std::stringstream tagname;
    tagname<<"Tag_";
    tagname<<tokens[0];
    TTree * sig = acceptedSig.tree(signame.str().c_str());
    sig->Write(signame.str().c_str());

 TTree * tagTree = acceptedTag.tree(tagname.str().c_str());
 tagTree->Write(tagname.str().c_str());

   // acceptedSig.tree(tokens[0])->Write("Signal");
 //  acceptedTag.tree(tokens[0])->Write("Tag");
  
  INFO( "Projecting" );
  
 std::vector<std::string> dalitzNames = {"01", "02", "12"}; 
  auto plots = acceptedSig.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  int i=0;
  for ( auto& plot : plots ) {
         std::stringstream projName;
         projName<<"Signal_vs_"<<tokens[0]<<"_s"<<dalitzNames[i];

         plot->Write(projName.str().c_str());
         i++;
    }
          auto proj = eventType.defaultProjections(nBins);      
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
          std::stringstream projName;
          projName<<"Signal_vs_"<<tokens[0]<<"_s"<<dalitzNames[i]<<"_vs_"<<dalitzNames[j];
        acceptedSig.makeProjection( Projection2D(proj[i], proj[j]), LineColor(kBlack) )->Write(projName.str().c_str()); 
      }
    } 

  if (outputVals){
    std::ofstream out;
    auto ampFile_tag = tokens[0] + "_" + ampFile;
    out.open(ampFile_tag.c_str());
    for (size_t i=0; i < acceptedSig.size(); i++){
      auto eventSig = acceptedSig[i];
      auto eventTag = acceptedTag[i];
      auto val = cs.getVals(eventSig, eventTag);
//      out<<ABCD<<"\n";
      auto sig01 = eventSig.s(0,1);
      auto sig02 = eventSig.s(0,2);
      auto sig12 = eventSig.s(1,2);
      auto tag01 = eventTag.s(0,1);
      auto tag02 = eventTag.s(0,2);
      auto tag12 = eventTag.s(1,2);
      auto A = val[0];
      auto B = val[1];
      auto C = val[2];
      auto D = val[3];
      auto ABCD = val[4];
      auto corr = val[5];
      out<<sig01<<"\t"<<sig02<<"\t"<<sig12<<"\t"
         <<A.real()<<"\t"<<A.imag()<<"\t"
         <<C.real()<<"\t"<<C.imag()<<"\t"
         <<tag01<<"\t"<<tag02<<"\t"<<tag12<<"\t"
         <<B.real()<<"\t"<<B.imag()<<"\t"
         <<D.real()<<"\t"<<D.imag()<<"\t"
         <<ABCD.real()<<"\t"<<ABCD.imag()<<"\t"
         <<corr.real()<<"\t"<<corr.imag()<<"\n";
         
    }
    out.close();
  }
INFO("End of Tag Loop");
 }
    f->Close();

 

 return 0;
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





EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string intFile){
   
    INFO("Getting Tag name split "<<tagName);
    auto tokens       = split(tagName, ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto types = makeEventTypes(sigName, tagName);
    INFO("Build EventType");
    auto signalType = types["signal"];
    auto tagType = types["tag"];
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");


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

    sigMCEvents = EventList(intFile + signame.str() , signalType);
    tagMCEvents = EventList(intFile + tagname.str() , tagType);
   
    if (type=="sigMC"){
        return sigMCEvents;
    }
    else if (type=="tagMC") {
        return tagMCEvents;
    }
    else {
        return null;
    }
    

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
