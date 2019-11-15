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
                       , TRandom* rndm )
{
  QcGenerator<PRIOR_TYPE> signalGenerator( priorSig, priorTag );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, eventsSig, eventsTag, nEvents );
}

void add_CP_conjugate( MinuitParameterSet& mps );

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 100, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
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



 for( auto& tag : tags ){


    auto tokens       = split(tag, ' ');
    INFO("tag = "<<tokens[0]);
 
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType();
  TRandom3 rand;
  rand.SetSeed( seed + 934534 );

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  INFO("makeCPConj = "<<makeCPConj);
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(MPS);
  }

  Particle p;
  bool debug=false;
  

//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );






//  EventType eventType2( NamedParameter<std::string>( "EventType2" , "", "EventType to generate second lot of events, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
//                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );



  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  INFO("Generating time-dependence? " << eventType.isTimeDependent() );
  EventList acceptedSig( sigType );
  EventList acceptedTag( tagType );

  INFO("Generating events with type = " << sigType << " and "<<tagType);


 // if ( genType == generatorType::CoherentSum ) 
    CorrelatedSum cs( sigType, tagType, MPS );
    PhaseSpace phspSig(sigType,&rand);
    PhaseSpace phspTag(tagType,&rand);
    INFO("Generating Events now!");
    GenerateEvents( acceptedSig, acceptedTag, cs, phspSig, phspTag , nEvents, blockSize, &rand );



      auto headPhsp = EventType({"psi(3770)0","D0","Dbar0"});
      auto psi_q = PhaseSpace(headPhsp); 
      auto beta  = [](const Event& event, const size_t&j){ return sqrt( event[4*j+0]*event[4*j+0] + event[4*j+1]*event[4*j+1] + event[4*j+2]*event[4*j+2] )/event[4*j+3] ; };
      auto p4     = [](const Event& event, const size_t&j){ return std::make_tuple(event[4*j+0], event[4*j+1], event[4*j+2]); };
    for (size_t i=0 ; i < acceptedSig.size(); i++){
        auto psi_event = psi_q.makeEvent();
        boost( acceptedSig[i], p4(psi_event,0), beta(psi_event,0));
        boost( acceptedTag[i]   , p4(psi_event,1), beta(psi_event,1));
      }


    if( acceptedSig.size() == 0 ) return -1;
    TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
    INFO( "Writing output file " );
    TTree * sig = acceptedSig.tree(tokens[0]);
   
    sig->Write("Signal");

 TTree * tagTree = acceptedTag.tree(tokens[0]);
 tagTree->Write("Tag");

   // acceptedSig.tree(tokens[0])->Write("Signal");
 //  acceptedTag.tree(tokens[0])->Write("Tag");
  
  
  
  auto plots = acceptedSig.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
    auto proj = eventType.defaultProjections(nBins);
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
        acceptedSig.makeProjection( Projection2D(proj[i], proj[j]), LineColor(kBlack) )->Write(); 
      }
    } 
    f->Close();
 }
 return 0;
}

/*
  DTEvent tmp(at(0).signal, at(0).tag);
  std::vector<int> id_sig(m_sigType.size()), 
    ids_sig(m_sigType.size()), 
    id_tag(m_tagType.size()),
    ids_tag(m_tagType.size());
 TTree* outputTree = new TTree(name.c_str(),name.c_str());
  for(size_t i = 0 ; i < sigType.size(); ++i )
  {
    outputTree->Branch((particleName(sigType, i)+"_PX").c_str(), &tmp.signal[4*i+0]); 
    outputTree->Branch((particleName(sigType, i)+"_PY").c_str(), &tmp.signal[4*i+1]); 
    outputTree->Branch((particleName(sigType, i)+"_PZ").c_str(), &tmp.signal[4*i+2]); 
    outputTree->Branch((particleName(sigType, i)+"_E").c_str(),  &tmp.signal[4*i+3]);
    outputTree->Branch((particleName(sigType, i)+"_ID").c_str(), &id_sig[i]);
    ids_sig[i] = ParticlePropertiesList::get( sigType[i] )->pdgID();
  }  
  for(size_t i = 0 ; i < m_tagType.size(); ++i )
  {
    outputTree->Branch(("Tag_"+particleName(tagType, i)+"_PX").c_str(), &tmp.tag[4*i+0]); 
    outputTree->Branch(("Tag_"+particleName(tagType, i)+"_PY").c_str(), &tmp.tag[4*i+1]); 
    outputTree->Branch(("Tag_"+particleName(tagType, i)+"_PZ").c_str(), &tmp.tag[4*i+2]); 
    outputTree->Branch(("Tag_"+particleName(tagType, i)+"_E").c_str(),  &tmp.tag[4*i+3]); 
    outputTree->Branch(("Tag_"+particleName(tagType, i)+"_ID").c_str(), &id_tag[i]);
    ids_tag[i] = ParticlePropertiesList::get( tagType[i] )->pdgID();
  }
  for( auto& evt: *this ){
    bool swap = gRandom->Uniform() > 0.5;
    tmp.set(evt.signal, evt.tag);
    if( swap ) tmp.invertParity();
    for(size_t i=0; i != m_sigType.size(); ++i)
      id_sig[i] = swap ? -ids_sig[i] : ids_sig[i];
    for(size_t i=0; i != m_tagType.size(); ++i)
      id_tag[i] = swap ? -ids_tag[i] : ids_tag[i];
    outputTree->Fill();
  }
    outputTree->Write();
    f->Close();
*/
  //} 
 // else {
 //   FATAL("Did not recognise configuration: " << genType );
 // }
 

 
  /*
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  acceptedSig.tree( "DalitzEventList" )->Write();
  auto plots = acceptedSig.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(nBins);
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){ 
        acceptedSig.makeProjection( Projection2D(proj[i], proj[j]), LineColor(kBlack) )->Write(); 
      }
    }
  } 
  */

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
        sgn = reOrIm == "Re" ? p.quasiCP() : 1; 
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





