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
#include "TEventList.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

using namespace AmpGen;

int main( int argc, char* argv[] )
{
  OptionsParser::setArgs( argc, argv );

  std::string inputFilename  = NamedParameter<std::string>( "Input", ( std::string ) "" ).getVal();
  std::string outputFilename = NamedParameter<std::string>( "Output", ( std::string ) "" ).getVal();
  std::string pdfLibrary     = NamedParameter<std::string>( "PdfLibrary", ( std::string ) "" ).getVal();
  std::string motherID       = NamedParameter<std::string>( "MotherID", ( std::string ) "" ).getVal();
  std::string treeName       = NamedParameter<std::string>( "Tree", ( std::string ) "" );
  std::string cuts           = "";
  unsigned int cIndex        = 0;
  std::string cut            = "";
  INFO( "Reading file " << inputFilename );

  TFile* f = TFile::Open( inputFilename.c_str(), "READ" );
  INFO( "Reading tree " << treeName );

  TTree* in_tree = (TTree*)f->Get( treeName.c_str() );
  in_tree->SetBranchStatus( "*", 1 );
  do {
    cut               = NamedParameter<std::string>( "Cuts" + std::to_string( cIndex++ ), ( std::string ) "" );
    double efficiency = in_tree->Draw( "", cut.c_str() ) / (double)in_tree->GetEntries();
    INFO( cut << " efficiency = " << efficiency );
    if ( cut != "" ) cuts += ( cuts == "" ? "" : "&&" ) + cut;
  } while ( cut != "" );

  INFO( "Using cut = " << cuts );
  bool usePIDCalib              = NamedParameter<int>( "usePIDCalib", 0 ).getVal();
  bool isMC                     = NamedParameter<int>( "isMC", 0 ).getVal();
  bool rejectMultipleCandidates = NamedParameter<int>( "rejectMultipleCandidates", 1 ).getVal();
  std::string plotsName         = NamedParameter<std::string>( "Plots", ( std::string ) "plots.root" ).getVal();

  if ( inputFilename == "" ) {
    ERROR( "No input specified in options" );
    return -1;
  }
  if ( treeName == "" ) {
    ERROR( "No tree specified in options" );
    return -1;
  }
  if ( outputFilename == "" ) {
    ERROR( "No output specified in options" );
    return -1;
  }
  if ( f == nullptr ) ERROR( inputFilename + " not found" );

  if ( in_tree == nullptr ) {
    ERROR( treeName + " not found" );
    return -1;
  } else
    INFO( "Got tree " << inputFilename << ":" << treeName );
  INFO( "Is MC ? " << ( isMC ? " True " : " False " ) );

  in_tree->Draw( ">>elist", cuts.c_str() );
  TEventList* elist = (TEventList*)gDirectory->Get( "elist" );
  INFO( "Total efficiency = " << elist->GetN() / (double)in_tree->GetEntries() );

  std::vector<size_t> eventsToTake;
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
        auto evtId = std::make_pair( eventNumber, runNumber );
        multipleCandidateIds[evtId].push_back( entry );
      }
    }

    std::mt19937 rng( 7 ); // random-number engine used (Mersenne-Twister in this case)
    // unsigned int multipleCandidatesRejected = 0;
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

  EventType evtType( NamedParameter<std::string>( "EventType" ).getVector() );
  

  std::vector<std::string> particles       = NamedParameter<std::string>( "ParticleNames" ).getVector();
  std::vector<std::string> monitorBranches = NamedParameter<std::string>( "Monitors" ).getVector();
  std::vector<std::string> branchFormat    = NamedParameter<std::string>( "BranchFormat" ).getVector();

  std::vector<std::string> branches;
  for ( auto& particle : particles ) {
    char buffer[100];
    sprintf( buffer, branchFormat[0].c_str(), particle.c_str() );
    branches.push_back( std::string( buffer ) );
    
    sprintf( buffer, branchFormat[1].c_str(), particle.c_str() );
    branches.push_back( std::string( buffer ) );
    
    sprintf( buffer, branchFormat[2].c_str(), particle.c_str() );
    branches.push_back( std::string( buffer ) );
    
    sprintf( buffer, branchFormat[3].c_str(), particle.c_str() );
    branches.push_back( std::string( buffer ) );
  }

  for ( auto& branch : monitorBranches ) branches.push_back( branch );

  INFO( "Branches = ["<< vectorToString(branches, ", " ) << "]" );

  in_tree->SetBranchStatus( "*", 0 );
  INFO("Constructing eventList");
  EventList evts( in_tree, evtType, Branches(branches), EntryList(eventsToTake) );

  if ( motherID != "" ) {
    INFO( "Converting D0bar -> D0" );
    in_tree->SetBranchStatus( "*", 0 );
    int id = 0;
    in_tree->SetBranchStatus( motherID.c_str() );
    in_tree->SetBranchAddress( motherID.c_str(), &id );
    for ( unsigned int i = 0; i < eventsToTake.size(); ++i ) {
      in_tree->GetEntry( eventsToTake[i] );
      if ( id < 0 ) {
        evts[i].invertParity( 4 );
      }
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
    };
  };

  TFile* outputFile = TFile::Open( outputFilename.c_str(), "RECREATE" );
  TTree* outputTree = evts.tree( "DalitzEventList");
  outputTree->Write();
  outputFile->Close();
  TFile* outputPlotFile = TFile::Open( plotsName.c_str(), "RECREATE" );
  auto plots            = evts.makeDefaultPlots();
  for ( auto& plot : plots ) {
    INFO( "Writing plot " << plot->GetName() << " to file" );
    plot->Write();
  }
  outputPlotFile->Close();
}
