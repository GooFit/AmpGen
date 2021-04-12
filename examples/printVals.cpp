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
#include "AmpGen/AddCPConjugate.h"
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


  
EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile){
   
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


    
    if (type=="signal"){
     return sigEvents;
    }
    else if (type=="tag"){ 
        return tagEvents;
    }
    else {
        return null;
    }
    


}


void printOut(EventList eventsSig, EventList eventsTag, pCorrelatedSum cs_tag, std::string output){
  std::ofstream out;
  out.open(output.c_str()); 
  int nElements = 6;

  TMatrixD * cov = new TMatrixD(nElements, nElements);
  std::string covFile = NamedParameter<std::string>("CovFile", "");

  if (covFile==""){
    for (int i=0;i<nElements;i++){
      for (int j =0;j<nElements;j++){
        (*cov)[i][j] = 0;
      }
    }
  }
 
  else{
    TFile * file = TFile::Open(covFile.c_str());
    TMatrixD * V = (TMatrixD*)file->Get("CovMatrix");
    cov = V;
    if (NamedParameter<bool>("pCorrelatedSum::printCov", false)){
    V->Print();
    }
    file->Close();
  }


  int printFreq = NamedParameter<int>("printFreq", 1000, "Print every n events");
  for (int i=0; i < eventsSig.size(); i++){
    //double percentage = (double)i/(double)eventsSig.size();
    //if (percentage != 0 && std::fmod(percentage, 0.01)==0)
    if(i%printFreq==0)INFO ("At "<<i<<"/"<<eventsSig.size());



    auto evt_sig = eventsSig[i];
    auto evt_tag = eventsTag[i];

    auto s01 = evt_sig.s(0,1);
    auto s02 = evt_sig.s(0,2);
    auto s12 = evt_sig.s(1,2);
    auto vals = cs_tag.getVals(evt_sig, evt_tag);

    auto aA = cs_tag.getA();
    auto aB = cs_tag.getB();
    auto aC = cs_tag.getC();
    auto aD = cs_tag.getD();

    auto A = aA.getVal(evt_sig);
    auto B = aB.getVal(evt_tag);
    auto C = aC.getVal(evt_sig);
    auto D = aD.getVal(evt_tag);

    auto A2 = aA(evt_sig)/eventsSig.size();
    auto B2 = aB(evt_tag)/eventsSig.size();
    auto C2 = aC(evt_sig)/eventsSig.size();
    auto D2 = aD(evt_tag)/eventsSig.size();


    auto ABCD = cs_tag.getVal(evt_sig, evt_tag);
    auto corr = vals[5];
    auto ACst = A * std::conj(C);
    auto dCorr = std::complex<real_t>(0,0);//cs_tag.errcorrection(evt_sig, cov);

    auto prob = cs_tag.prob(evt_sig, evt_tag)/eventsSig.size();
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
}



void printPerBin(int nBins, pCorrelatedSum cs_tag, EventList sigevents_tag, EventList tagevents_tag, std::string tagName){
  int nElements = 6;

  TMatrixD * cov = new TMatrixD(nElements, nElements);
  std::string covFile = NamedParameter<std::string>("CovFile", "");

  if (covFile==""){
    for (int i=0;i<nElements;i++){
      for (int j =0;j<nElements;j++){
        (*cov)[i][j] = 0;
      }
    }
  }
 
  else{
    TFile * file = TFile::Open(covFile.c_str());
    TMatrixD * V = (TMatrixD*)file->Get("CovMatrix");
    cov = V;
    if (NamedParameter<bool>("pCorrelatedSum::printCov", false)){
    V->Print();
    }
    file->Close();
  }



for (int j=1;j<nBins + 1;j++){
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
    auto ACst = A * std::conj(C);
    auto dCorr =  std::complex<real_t>(0,0);//cs_tag.errcorrection(evt_sig, cov);
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


void printDT(int i, int j, EventList sigEvents, EventList tagEvents, pCorrelatedSum cs_tag, std::string pref=""){
  std::ofstream out;
  std::stringstream outName;
  outName << "Kspipi"<<pref<<i<<j<<".csv";
  out.open(outName.str().c_str());


  std::ofstream outPM;
  std::stringstream outNamePM;
  outNamePM << "Kspipi"<<pref<<i<<-j<<".csv";
  outPM.open(outNamePM.str().c_str());


  std::ofstream outMP;
  std::stringstream outNameMP;
  outNameMP << "Kspipi"<<pref<<-i<<j<<".csv";
  outMP.open(outNameMP.str().c_str());
 
 

  std::ofstream outMM;
  std::stringstream outNameMM;
  outNameMM << "Kspipi"<<pref<<-i<<-j<<".csv";
  outMM.open(outNameMM.str().c_str());
  int nElements = 6;

  TMatrixD * cov = new TMatrixD(nElements, nElements);
  std::string covFile = NamedParameter<std::string>("CovFile", "");

  if (covFile==""){
    for (int i=0;i<nElements;i++){
      for (int j =0;j<nElements;j++){
        (*cov)[i][j] = 0;
      }
    }
  }
 
  else{
    TFile * file = TFile::Open(covFile.c_str());
    TMatrixD * V = (TMatrixD*)file->Get("CovMatrix");
    cov = V;
    if (NamedParameter<bool>("pCorrelatedSum::printCov", false)){
    V->Print();
    }
    file->Close();
  }


  

  auto dd_Low_i = 2 * M_PI * i/8. - 5 * M_PI/4.;
  auto dd_Low_j = 2 * M_PI * j/8. - 5 * M_PI/4.;
  auto dd_High_i = 2 * M_PI * (i+1)/8. - 5 * M_PI/4.;
  auto dd_High_j = 2 * M_PI * (j+1)/8. - 5 * M_PI/4.;

  for (int evtNum=0; evtNum<sigEvents.size(); evtNum++){
    auto evt_sig = sigEvents[evtNum];
    auto evt_tag = tagEvents[evtNum];
    auto sig01 = evt_sig.s(0,1);
    auto sig02 = evt_sig.s(0,2);
    auto sig12 = evt_sig.s(1,2);

    auto tag01 = evt_tag.s(0,1);
    auto tag02 = evt_tag.s(0,2);
    auto tag12 = evt_tag.s(1,2);
 

    auto vals = cs_tag.getVals(evt_sig, evt_tag);

    auto aA = cs_tag.getA();
    auto aB = cs_tag.getB();
    auto aC = cs_tag.getC();
    auto aD = cs_tag.getD();

    auto A = aA.getVal(evt_sig);
    auto B = aB.getVal(evt_tag);
    auto C = aC.getVal(evt_sig);
    auto D = aD.getVal(evt_tag);


    auto ABCD = cs_tag.getVal(evt_sig, evt_tag);
    auto corr = vals[5];
    auto prob = cs_tag.prob(evt_sig, evt_tag);

    auto ACst = A * std::conj(C);
    auto DBst = D * std::conj(B);
    auto dCorr =  std::complex<real_t>(0,0);//cs_tag.errcorrection(evt_sig, cov);
    auto A2 = aA(evt_sig);
    auto B2 = aB(evt_tag);
    auto C2 = aC(evt_sig);
    auto D2 = aD(evt_tag);



    auto cosDD = -((fcn::pow(abs(ABCD), 2)() - fcn::pow(abs(A), 2)() * fcn::pow(abs(B), 2)() - fcn::pow(abs(C), 2)() * fcn::pow(abs(D), 2)())/(2 * abs(A) * abs(B) * abs(C) * abs(D) )).real();

    auto dd_sig = std::arg(ACst);
    auto dd_tag = std::arg(DBst);
    auto ACstCorr = A * std::conj(C) * fcn::exp(Constant(0, 1) * corr)();
    auto ddCorr = std::arg(ACstCorr);
 
    bool cond_sig = (dd_sig > dd_Low_i) && (dd_sig < dd_High_i);
    bool cond_tag = (dd_tag > dd_Low_j) && (dd_tag < dd_High_j);

    if (cond_sig && cond_tag){
      bool PP = (sig01 > sig02) && (tag01>tag02);
      bool PM = (sig01 > sig02) && (tag01<tag02);
      bool MP = (sig01 < sig02) && (tag01>tag02);
      bool MM = (sig01 < sig02) && (tag01<tag02);
      if (PP) out << sig01 << "\t" << sig02 << "\t"<< tag01 << "\t" <<tag02 << "\t" << dd_sig << "\t"<<dd_tag<<"\t"<<ABCD.real()<<"\t"<<ABCD.imag()<<"\t"<<A.real()<<"\t"<<A.imag()<<"\t"<<B.real()<<"\t"<<B.imag()<<"\t"<<C.real()<<"\t"<<C.imag()<<"\t"<<D.real()<<"\t"<<D.imag()<<"\n";
      if (PM) outPM << sig01 << "\t" << sig02 << "\t"<< tag01 << "\t" <<tag02 << "\t" << dd_sig << "\t"<<dd_tag<<"\t"<<ABCD.real()<<"\t"<<ABCD.imag()<<"\t"<<A.real()<<"\t"<<A.imag()<<"\t"<<B.real()<<"\t"<<B.imag()<<"\t"<<C.real()<<"\t"<<C.imag()<<"\t"<<D.real()<<"\t"<<D.imag()<<"\n";
      if (MP) outMP << sig01 << "\t" << sig02 << "\t"<< tag01 << "\t" <<tag02 << "\t" << dd_sig << "\t"<<dd_tag<<"\t"<<ABCD.real()<<"\t"<<ABCD.imag()<<"\t"<<A.real()<<"\t"<<A.imag()<<"\t"<<B.real()<<"\t"<<B.imag()<<"\t"<<C.real()<<"\t"<<C.imag()<<"\t"<<D.real()<<"\t"<<D.imag()<<"\n";
      if (MM) outMM << sig01 << "\t" << sig02 << "\t"<< tag01 << "\t" <<tag02 << "\t" << dd_sig << "\t"<<dd_tag<<"\t"<<ABCD.real()<<"\t"<<ABCD.imag()<<"\t"<<A.real()<<"\t"<<A.imag()<<"\t"<<B.real()<<"\t"<<B.imag()<<"\t"<<C.real()<<"\t"<<C.imag()<<"\t"<<D.real()<<"\t"<<D.imag()<<"\n";
    }
  }
  out.close();
  outPM.close();
  outMP.close();
  outMM.close();
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
  bool doQC  = NamedParameter<bool>("doQC", true, "Use QC");
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
  int NInt = NamedParameter<int>("NInt", 10000000);
  int NPoints = NamedParameter<int>("NPoints", 10000);


  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile = NamedParameter<std::string>("IntegrationSample", ""          , "Name of file containing flat sample." );




  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

if (doQC){
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    //add_CP_conjugate(MPS);
    AddCPConjugate(MPS);
  }
  for (int i=0; i < tags.size(); i++){



    INFO("Getting Tag name split "<<tags[i]);
    auto tokens       = split(tags[i], ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto types = makeEventTypes(pNames, tags[i]);
    INFO("Build EventType");
    auto signalType = types["signal"];
    auto tagType = types["tag"];


    EventList sigevents_tag = Generator<>(signalType, &rndm).generate(NPoints);
    EventList tagevents_tag = Generator<>(tagType, &rndm).generate(NPoints);


    if (dataFile != ""){
      sigevents_tag = getEvents("signal", pNames, tags[i], dataFile);
      tagevents_tag = getEvents("tag", pNames, tags[i], dataFile);
    }


 


    EventList sigMCevents_tag = Generator<>(sigType, &rndm).generate(NInt);
    EventList tagMCevents_tag = Generator<>(tagType, &rndm).generate(NInt);
 
    CoherentSum sig(signalType, MPS);
    pCorrelatedSum cs_tag (signalType, tagType, MPS);

     
    cs_tag.setEvents(sigevents_tag, tagevents_tag);
    cs_tag.setMC(sigMCevents_tag, tagMCevents_tag);
    cs_tag.prepare();
    std::ofstream out; 
    
    sig.setEvents(sigevents_tag);
    sig.setMC(sigMCevents_tag);
    sig.prepare();
    auto evt = sigevents_tag[0];
    auto evtT = evt; evtT.swap(1,2);
    auto v = sig.getValNoCache(evt);
    auto vT =  sig.getValNoCache(evtT);
    INFO("psi("<<evt.s(0,1)<<", "<<evt.s(0,2)<<") = "<<v);
    INFO("psi("<<evtT.s(0,1)<<", "<<evtT.s(0,2)<<") = "<<vT);


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
    std::stringstream tag_logMC;
    tag_logMC<<tagName<<"MC_vals.csv";

    std::ofstream outMC;
    outMC.open(tag_logMC.str().c_str()); 

    printOut(sigevents_tag, tagevents_tag, cs_tag, tag_log.str());
    //printOut(sigMCevents_tag, tagMCevents_tag, cs_tag, tag_logMC.str());

    std::stringstream tagNameMC;
    tagNameMC<<tagName<<"MC";

    if (NamedParameter<bool>("BinPrint", false)){ 
    printPerBin(8, cs_tag, sigevents_tag, tagevents_tag, tagName);
    //printPerBin(8, cs_tag, sigMCevents_tag, tagMCevents_tag, tagNameMC.str());
    if (tagName=="Kspipi"){
      for (int i=1;i<9;i++){
        for (int j=1;j<9;j++){
        printDT(i, j, sigevents_tag, tagevents_tag, cs_tag);
        //printDT(i,j, sigMCevents_tag, tagMCevents_tag, cs_tag, "MC");
        }
      }
    }
  }


  }
}

else{
   MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(MPS);
  }
   EventList events(dataFile, eventType);

  EventList eventsMC =  Generator<>(eventType, &rndm).generate(NInt);
  CoherentSum sig(eventType, MPS);
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
