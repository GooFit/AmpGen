#include <string>
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventType.h"
#include "AmpGen/EventList.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Particle.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/Generator.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/ArgumentPack.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

using namespace std;
using namespace AmpGen;

std::vector<std::string> makeBranches(EventType Type, std::string prefix);
void add_CP_conjugate( MinuitParameterSet& mps );


int main(int argc, char** argv ){
    OptionsParser::setArgs( argc, argv );

    INFO("Test for CorrelatedSum");
    auto pNames = NamedParameter<std::string>("EventType" , "D0 K0S0 pi+ pi-"    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
    auto tags = NamedParameter<std::string>("TagTypes" , std::string(),
     "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
    std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
    bool doCorrFit  = NamedParameter<bool>("doCorrFit", true, "Fit the correlated pdf");
    bool doDebugNorm  = NamedParameter<bool>("doDebugNorm", false, "Debug the normalisation of the pdf");
    bool doAFit  = NamedParameter<bool>("doAFit", false, "Fit the A pdf in AB-CD correlated pdf");
    bool doBFit  = NamedParameter<bool>("doBFit", false, "Fit the B pdf in BB-CD correlated pdf");
    bool doCFit  = NamedParameter<bool>("doCFit", false, "Fit the C pdf in CB-CD correlated pdf");
    bool doDFit  = NamedParameter<bool>("doDFit", false, "Fit the D pdf in DB-CD correlated pdf");
    bool m_debug = NamedParameter<bool>("debug", false, "Debug flag");
    MinuitParameterSet mps;    
    mps.loadFromStream();
    //add_CP_conjugate(mps);
    TFile * data = TFile::Open(dataFile.c_str());
    TFile * mc = TFile::Open(intFile.c_str());
    TRandom3 rndm;
    int seed=0;
    rndm.SetSeed( seed );
    gRandom = &rndm;
    EventType signalType(pNames);
    EventList sigEvents(signalType);
    EventList sigMCEvents(signalType);
//    INFO("The EventType is "<<pNames); 
    EventList sigonlyEvents = Generator<>(signalType, &rndm).generate(2e5) ;
    EventList sigonlyMCEvents = Generator<>(signalType, &rndm).generate(2e6) ;


    for (auto tag : tags){
        INFO("Tag = "<<tag );
        auto tokens       = split(tag, ' ');
        auto tagParticle  = Particle(tokens[1], {}, false);
        EventType tagType = tagParticle.eventType();
        EventList tagonlyEvents = Generator<>(tagType, &rndm).generate(2e5) ;
        EventList tagonlyMCEvents = Generator<>(tagType, &rndm).generate(2e6) ;
        EventList tagEvents(tagType);
        EventList tagMCEvents(tagType);

        TTree * dataTree = (TTree*)data->Get(tokens[0].c_str());
        TTree * mcTree = (TTree*)mc->Get(tokens[0].c_str());
        auto sigBranches = makeBranches(signalType, "");
        auto tagBranches = makeBranches(tagType, "Tag_");
        auto argPackSig = ArgumentPack(Branches(sigBranches));
        auto argPackTag = ArgumentPack(Branches(tagBranches));

           
        INFO("Loading Tag Events from DataSample");
        tagEvents.loadFromTree(dataTree, argPackTag);
        INFO("Loading Tag Events from IntegrationSample");
        tagMCEvents.loadFromTree(mcTree, argPackTag);

        INFO("Loading Signal Events from DataSample");
        sigEvents.loadFromTree(dataTree, argPackSig);

        INFO("Loading Signal Events from IntegrationSample");
        sigMCEvents.loadFromTree(mcTree, argPackSig);


        if (doAFit){

          CoherentSum A(signalType, mps);
          A.setEvents(sigonlyEvents);
          A.setMC(sigonlyMCEvents);
          A.prepare();
          INFO("norm A = "<<A.norm());
          auto ALL = make_pdf(A);
          ALL.setEvents(sigEvents);
          Minimiser miniA(ALL, &mps);
          auto fitA = miniA.doFit();
          //INFO("A Fit" << fitA);
          auto coVarA = miniA.covMatrix();
          INFO("Printing CoVariant matrix for A");
          coVarA.Print();
          INFO("LL(A) = "<<ALL.getVal());
        }

        if (doBFit){

          CoherentSum B(tagType.conj(true), mps);
          B.setEvents(tagonlyEvents);
          B.setMC(tagonlyMCEvents);
          B.prepare();
          INFO("norm B = "<<B.norm());
          auto BLL = make_pdf(B);
          BLL.setEvents(tagEvents);
          Minimiser miniB(BLL, &mps);
          auto fitB = miniB.doFit();
          INFO("B Fit" << fitB);
          INFO("LL(B) = "<<BLL.getVal());
        }

        if (doCFit){

          CoherentSum C(signalType.conj(true), mps);
          C.setEvents(sigonlyEvents);
          C.setMC(sigonlyMCEvents);
          C.prepare();
          INFO("norm C = "<<C.norm());
          auto CLL = make_pdf(C);
          CLL.setEvents(sigEvents);
          Minimiser miniC(CLL, &mps);
          auto fitC = miniC.doFit();
          INFO("C Fit" << fitC);
          INFO("LL(C) = "<<CLL.getVal());
        }

        if (doDFit){

          CoherentSum D(tagType, mps);
          D.setEvents(tagonlyEvents);
          D.setMC(tagonlyMCEvents);
          D.prepare();
          INFO("norm D = "<<D.norm());
          auto DLL = make_pdf(D);
          DLL.setEvents(tagEvents);
          Minimiser miniD(DLL, &mps);
          auto fitD = miniD.doFit();
          INFO("D Fit" << fitD);
          INFO("LL(D) = "<<DLL.getVal());
        }


        CorrelatedSum cs(signalType, tagType, mps);
        cs.setEvents(sigEvents, tagEvents);
        cs.setMC(sigMCEvents, tagMCEvents);
        cs.prepare();
       // if (doDebugNorm) cs.debugNorm();
        
        
        if (doCorrFit){
          auto i=0;
          //N fits - idea is to see how "stable" our fit is - a simple test
          //just redoes the fit n times, obviously should be identical results!
          auto n=1;
          while(i < n){
            auto csLL = make_likelihood(sigEvents, tagEvents,  cs);
  //          csLL.setEvents(sigEvents, tagEvents);
            //auto LL = csLL.getVal();
            INFO("Beginning the fit"); 
            Minimiser mini(csLL, &mps);
            auto fcn = mini.FCN();
            INFO("fcn = "<<fcn);
            //INFO("LL = "<<LL);
            mini.prepare();
        
            if (m_debug) mini.gradientTest();
          

            mini.doFit();
            auto covar = mini.covMatrix();
            INFO("Printing CoVariant Matrix");
            covar.Print();
            //auto covarFull = mini.covMatrixFull();
            //INFO("Printing Full CoVariant Matrix");
            //covarFull.Print();
            for_each(csLL.pdfs(), [&] (auto& pdf){
           //   corrEventList cEL_data(sigEvents, tagEvents);
           //   corrEventList cEL_mc(sigMCEvents, tagMCEvents);
            //  std::vector<std::vector<TH1D*> > hists_mc =  cEL_mc.makeDefaultProjections(WeightFunction(pdf), "test");
            });
            i++;
            }


          //INFO("Has the fit been done? "<<fit);
          INFO("Resetting the correlated sum");
        }
        cs.reset(true);
        sigEvents.resetCache();
        sigMCEvents.resetCache();
        tagEvents.resetCache();
        tagMCEvents.resetCache();



    }
    return 0;

}


std::vector<std::string> makeBranches(EventType Type, std::string prefix){
  auto n = Type.finalStates().size();
  std::vector<std::string> branches;
  std::vector<std::string> varNames = {"E", "PX", "PY", "PZ"};
  for (long unsigned int i=0; i<n; i++){
    auto part = replaceAll(Type.finalStates()[i], "+", "p");
    part = replaceAll(part, "-", "m");
    for (auto j: varNames){
      std::ostringstream stringStream;
      stringStream<<prefix<<part<<"_"<<j;
      branches.push_back(stringStream.str());
    }
  }
  return branches;
}

void add_CP_conjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  for( auto& param : mps ){
    const std::string name = param->name();
    size_t pos=0;
    std::string new_name = name; 
    int sgn=1;
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
      tmp.push_back( new MinuitParameter(new_name, Flag::Free, sgn * param->mean(), param->err(), 0, 0));
    }
  }
  for( auto& p : tmp ) mps.add( p );
}
