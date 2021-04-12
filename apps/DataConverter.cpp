#include <RtypesCore.h>
#include <TDirectory.h>
#include <algorithm>
#include <dlfcn.h>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/DynamicFCN.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Projection.h"

#include "TEventList.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

using namespace AmpGen;

void invertParity( Event& event, const size_t& nParticles=0)
{
  for( size_t i = 0 ; i < nParticles; ++i )
  {
    event[4*i + 0] = -event[4*i+0];
    event[4*i + 1] = -event[4*i+1];
    event[4*i + 2] = -event[4*i+2];
  }
}


int main( int argc, char* argv[] )
{
  OptionsParser::setArgs( argc, argv );
  std::string inputFilename                = NamedParameter<std::string>("Input"     , "", "Input ROOT file" );
  std::string treeName                     = NamedParameter<std::string>("Tree"      , "", "Input ROOT tree." );
  std::string outputFilename               = NamedParameter<std::string>("Output"    , "", "Output ROOT file" );
  std::string pdfLibrary                   = NamedParameter<std::string>("PdfLibrary", "", "PDF Library that used to generate this sample for MC reweighting (MC only)" );
  std::string motherID                     = NamedParameter<std::string>("MotherIDBranch"    , "", "Name of branch that contains the ID of the parent, i.e. > 0 for particles, < 0 for antiparticles." );
  std::string plotsName                    = NamedParameter<std::string>("Plots"     ,"plots.root", "Output file for ROOT plots");
  std::vector<std::string> particles       = NamedParameter<std::string>("ParticleNames" , std::vector<std::string>() ).getVector();
  std::vector<std::string> monitorBranches = NamedParameter<std::string>("Monitors"      , std::vector<std::string>() ).getVector();
  std::vector<std::string> branchFormat    = NamedParameter<std::string>("BranchFormat"  , std::vector<std::string>() ).getVector();
  std::vector<std::string> friends         = NamedParameter<std::string>("Friends"       , std::vector<std::string>() ).getVector();
  bool usePIDCalib                         = NamedParameter<bool>("usePIDCalib"             , false );
  bool rejectMultipleCandidates            = NamedParameter<bool>("rejectMultipleCandidates", true  );
  std::string cuts                         = vectorToString( NamedParameter<std::string>("Cut","").getVector() , " && "); 
  EventType evtType( NamedParameter<std::string>( "EventType" ).getVector() );

  INFO( "Reading file " << inputFilename );
  INFO( "Outputting file: " << outputFilename);
  TFile* f = TFile::Open( inputFilename.c_str(), "READ" );
  INFO( "Reading tree " << treeName );

  TTree* in_tree = (TTree*)f->Get( treeName.c_str() );
  in_tree->SetBranchStatus( "*", 1 );
  for( auto& frie : friends ){
    auto tokens = split( frie, ':');
    in_tree->AddFriend( tokens[1].c_str(), tokens[0].c_str() );
  }

  INFO( "Using cut = " << cuts );

  if(inputFilename  == "") FATAL("No input specified in options" );
  if(treeName       == "") FATAL("No tree specified in options" );
  if(outputFilename == "") FATAL("No output specified in options" );
  if(f       == nullptr  ) FATAL(inputFilename + " not found" );
  if(in_tree == nullptr  ) FATAL(treeName + " not found" );
  
  INFO( "Got tree " << inputFilename << ":" << treeName );

  in_tree->Draw( ">>elist", cuts.c_str() );
  TEventList* elist = (TEventList*)gDirectory->Get( "elist" );
  INFO( "Total efficiency = " << elist->GetN() / (double)in_tree->GetEntries() );

  std::vector<size_t> eventsToTake;
  
  std::vector<std::string> branches;
  for( auto& particle : particles ) {
    for(size_t i = 0 ; i < 4; ++i ) branches.push_back( mysprintf(branchFormat[i], particle.c_str()));
  }
  for(auto& branch : monitorBranches) branches.push_back( branch );
  
  if ( rejectMultipleCandidates ) {
    ULong64_t totCandidate;
    ULong64_t eventNumber;
    UInt_t runNumber;
    std::map<std::pair<ULong64_t, UInt_t>, std::vector<unsigned int>> multipleCandidateIds;

    in_tree->SetBranchStatus( "*", 0 );
    in_tree->SetBranchStatus( "totCandidates", 1 );
    in_tree->SetBranchAddress( "totCandidates", &totCandidate );
    in_tree->SetBranchStatus( "eventNumber", 1 );
    in_tree->SetBranchAddress( "eventNumber", &eventNumber );
    in_tree->SetBranchStatus( "runNumber", 1 );
    in_tree->SetBranchAddress( "runNumber", &runNumber );

    INFO( "Doing multiple candidate rejection ...." );
    for ( int i = 0; i < elist->GetN(); ++i ) {
      if ( i % 100000 == 0 ) INFO( "Processed: " << i << " events" );
      unsigned int entry = elist->GetEntry( i );
      in_tree->GetEntry( entry );
      if ( totCandidate == 1 ) {
        eventsToTake.push_back( entry );
      } else {
        auto evtId = std::make_pair(eventNumber, runNumber);
        multipleCandidateIds[evtId].push_back( entry );
      }
    }

    std::mt19937 rng(7); // random-number engine used (Mersenne-Twister in this case)
    for ( auto& candidate : multipleCandidateIds ) {
      unsigned int nCand = candidate.second.size();
      unsigned int j     = nCand == 1 ? 0 : std::uniform_int_distribution<int>( 0, nCand - 1 )( rng );
      eventsToTake.push_back( candidate.second[j] );
    }
    std::sort( eventsToTake.begin(), eventsToTake.end() );
    INFO( "Events before multiple candidate rejection = " << elist->GetN() << " after = " << eventsToTake.size() );
  } else {
    for ( int i = 0; i < elist->GetN(); ++i ) {
      eventsToTake.push_back( elist->GetEntry( i ) );
    }
  }

  EventList evts( in_tree, evtType, Branches(branches), EntryList(eventsToTake), GetGenPdf(false), ApplySym(true) );

  INFO( "Branches = ["<< vectorToString(branches, ", " ) << "]" );

  in_tree->SetBranchStatus( "*", 0 );
  INFO("Constructing eventList");

  if ( motherID != "" ) {
    INFO( "Converting " << evtType.mother() << " " << eventsToTake.size() << " " << evts.size()  );
    in_tree->SetBranchStatus( "*", 0 );
    int id = 0;
    in_tree->SetBranchStatus( motherID.c_str() );
    in_tree->SetBranchAddress( motherID.c_str(), &id );
    for ( unsigned int i = 0; i < eventsToTake.size(); ++i ) {
      in_tree->GetEntry( eventsToTake[i] );
      if ( id < 0 ) invertParity( evts[i] , evtType.size() );
    }
  }

  INFO( "Done building event list, size = " << evts.size() );

  if ( usePIDCalib ) {
    INFO( "Getting event weights from PID calib" );
    std::string stub_path = inputFilename.substr( 0, inputFilename.find_last_of( '/' ) );
    INFO( stub_path );
    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    std::vector<Float_t> weights( particles.size(), 0 );
    for ( unsigned int i = 0; i < particles.size(); ++i ) {
      files.push_back( TFile::Open( ( stub_path + "/pidCalib_" + particles[i] + "_repacked.root" ).c_str() ) );
      trees.push_back( (TTree*)( *files.rbegin() )->Get( "CalibTool_PIDCalibTree" ) );
      ( *trees.rbegin() )->SetBranchAddress( ( particles[i] + "_PIDCalibEffWeight" ).c_str(), &( weights[i] ) );
    }
    for ( unsigned int i = 0; i < eventsToTake.size(); ++i ) {
      double weight    = 1;
      unsigned int evt = eventsToTake[i];
      for ( auto& t : trees ) {
        if ( evt > t->GetEntries() ) {
          ERROR( "Accessing out of bounds data - something has gone HORRIBLY wrong" );
        }
        t->GetEntry( evt );
      }
      for ( auto& w : weights ) {
        if ( w < 0 ) {
          ERROR( "PID weight for event " << evt << " = " << w );
          evts[i].print();
        }
        weight *= w;
      }
      evts[i].setWeight( weight );
    }
  }
  evts.transform( [=](auto& event){ for( size_t i = 0 ; i < 4*evtType.size(); ++i ) event[i] /= 1000. ; } );

  if ( pdfLibrary != "" ) {
    INFO( "Setting generator level PDF from " << pdfLibrary );
    void* handle = dlopen( pdfLibrary.c_str(), RTLD_NOW );
    if ( handle == nullptr ) dlerror();

    DynamicFCN<double( const double*, const int& )> fcn( handle, "FCN" );
    for ( unsigned int i = 0; i < evts.size(); ++i ) {
      if ( i % 500000 == 0 ) {
        INFO( "Set for " << i << " events" );
      }
      evts[i].setGenPdf( fcn( (const real_t*)(evts[i]), 1 ) );
    }
  }
  INFO("Writing file: " << outputFilename);
  TFile* outputFile = TFile::Open( outputFilename.c_str(), "RECREATE" );
  INFO("Made file :-> " );

  TTree* outputTree = evts.tree("DalitzEventList");
  outputTree->Write();
  
  INFO("Closing file...");
  outputFile->Close();
  TFile* outputPlotFile = TFile::Open( plotsName.c_str(), "RECREATE" );
  auto plots            = evts.makeDefaultProjections();
  for ( auto& plot : plots ) {
    INFO( "Writing plot " << plot->GetName() << " to file" );
    plot->Write();
  }
  outputPlotFile->Close();
}
