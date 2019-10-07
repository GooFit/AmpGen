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


std::vector<std::string> makeBranches(EventType Type, std::string prefix){
  auto n = Type.finalStates().size();
  std::vector<std::string> branches;
  std::vector<std::string> varNames = {"PX", "PY", "PZ", "E"};
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
    std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
    std::string dataProj = NamedParameter<std::string>("dataProj", "data.csv", "name of csv file for data projection");
    std::string mcProj = NamedParameter<std::string>("mcProj", "mc.csv", "name of csv file for mc projection");
    std::string sigOnlyProj = NamedParameter<std::string>("sigOnlyProj", "sigOnly.csv", "name of csv file for sigOnly projection");
    std::string logFile  = NamedParameter<std::string>("LogFile"   , "QCFitter.log", "Name of the output log file");
    int nBins = NamedParameter<int> ("nBins", 100, "number of bins for the projection");

    bool m_debug = NamedParameter<bool>("debug", false, "Debug flag");   
    MinuitParameterSet mps;    
    mps.loadFromStream();
    //add_CP_conjugate(mps);

    for (auto tag : tags){
        INFO("Tag = "<<tag );
        TFile * data = TFile::Open(dataFile.c_str());
        TFile * mc = TFile::Open(intFile.c_str());
        TRandom3 rndm;
        int seed=0;
        rndm.SetSeed( seed );
        gRandom = &rndm;
        EventType signalType(pNames);
        EventList sigEvents(signalType);
        EventList sigMCEvents(signalType);
        EventList sigonlyEvents = Generator<>(signalType, &rndm).generate(2e5) ;


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

        INFO("Making Projections");
        
        TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
        dataTree->Write("DataTree");
        sigEvents.tree("Data")->Write();
        auto proj1 = signalType.defaultProjections(nBins);
        
        auto plots = sigEvents.makeProjections( proj1, LineColor(kBlack) ); 
        auto plotsSO = sigonlyEvents.makeProjections( proj1, LineColor(kBlack) ); 
       
        
        /*
        ofstream outfile;
        stringstream stream;
        stream<<dataProj;



        auto string = stream.str();
        outfile.open(string.c_str());

        for (auto& evt : sigEvents){
            outfile<<evt.s(1,0)<<"\t"<<evt.s(2,0)<<"\t"<<evt.s(2,1)<<"\t"<<evt.s(0,0)<<"\t"<<evt.s(1,1)<<"\t"<<evt.s(2,2)<<"\t"<<evt.weight()/evt.genPdf()<<"\n";
        }
        outfile.close();


        ofstream outfileMC;
        stringstream streamMC;
        streamMC<<mcProj;

        auto stringMC = streamMC.str();
        outfileMC.open(stringMC.c_str());

        for (auto& evt : sigMCEvents){
            evt.print();
            outfileMC<<evt.s(1,0)<<"\t"<<evt.s(2,0)<<"\t"<<evt.s(2,1)<<"\t"<<evt.s(0,0)<<"\t"<<evt.s(1,1)<<"\t"<<evt.s(2,2)<<"\t"<<evt.weight()/evt.genPdf()<<"\n";

            
        }
        outfileMC.close();

        ofstream outfilesigOnly;
        stringstream streamsigOnly;
        streamsigOnly<<sigOnlyProj;

        auto stringsigOnly = streamsigOnly.str();
        outfilesigOnly.open(stringsigOnly.c_str());

        for (auto& evt : sigonlyEvents){

            outfilesigOnly<<evt.s(1,0)<<"\t"<<evt.s(2,0)<<"\t"<<evt.s(2,1)<<"\t"<<evt.s(0,0)<<"\t"<<evt.s(1,1)<<"\t"<<evt.s(2,2)<<"\t"<<evt.weight()/evt.genPdf()<<"\n";

            
        }
        outfilesigOnly.close();
    
      */
       
        
        mcTree->Write("MCTree");
        sigMCEvents.tree("MC")->Write();
        auto sigplots = sigEvents.makeDefaultProjections(Prefix("sig_data"), Bins(100));
        auto sigplotsMC = sigMCEvents.makeDefaultProjections(Prefix("sig_mc"), Bins(100));
        auto tagplots = tagEvents.makeDefaultProjections(Prefix("tag_data"), Bins(100));
        auto tagplotsMC = tagMCEvents.makeDefaultProjections(Prefix("tag_mc"), Bins(100));

           
        for ( auto& plot : sigplots ) plot->Write();
        for ( auto& plot : sigplotsMC ) plot->Write();
        for ( auto& plot : tagplots ) plot->Write();
        for ( auto& plot : tagplotsMC ) plot->Write();



        output->Close();
        data->Close();
        mc->Close();
           

    }
    return 0;
}