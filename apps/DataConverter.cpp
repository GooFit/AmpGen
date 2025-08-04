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
#include "AmpGen/Property.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Projection.h"
#include "AmpGen/TreeReader.h"
#include "TEventList.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

using namespace AmpGen;

void invertParity(Event &event, const size_t &nParticles = 0)
{
  for(size_t i = 0; i < nParticles; ++i)
    {
      event[4 * i + 0] = -event[4 * i + 0];
      event[4 * i + 1] = -event[4 * i + 1];
      event[4 * i + 2] = -event[4 * i + 2];
    }
}

int main(int argc, char *argv[])
{
  using strings = std::vector<std::string>; 
  OptionsParser::setArgs(argc, argv);
  Property<std::string> inputFilename   {nullptr, "Input" , "", "Input ROOT file(s)"};
  Property<std::string> treeName        {nullptr, "Tree"  , "", "Input ROOT tree."};
  Property<std::string> outputFilename  {nullptr, "Output", "", "Output ROOT file"};
  Property<std::string> pdfLibrary      {nullptr, "PdfLibrary", "", "PDF Library that used to generate this sample for MC reweighting (MC only)"};
  Property<std::string> motherID        {nullptr, "MotherIDBranch", "", "Name of branch that contains the ID of the parent, i.e. > 0 for particles, < 0 for antiparticles."};
  Property<std::string> plotsName       {nullptr, "Plots", "plots.root", "Output file for ROOT plots"};
  Property<strings> particles           {nullptr, "ParticleNames"}; 
  Property<strings> monitorBranches     {nullptr, "Monitors"}; 
  Property<strings> branchFormat        {nullptr, "BranchFormat"}; 
  Property<strings> friends             {nullptr, "Friends"}; 
  Property<strings> idBranches          {nullptr, "IdBranches", std::vector<std::string>()}; 
  Property<std::string> units           {nullptr, "Units", "MeV"};
  Property<bool> usePIDCalib            {nullptr, "usePIDCalib", false};
  Property<bool> rejectMultipleCandidates {nullptr, "rejectMultipleCandidates", true};
  Property<strings> cuts                {nullptr, "Cut"};
  Property<strings> evtType_s           {nullptr, "EventType"}; 
  EventType evtType(evtType_s);

  std::vector<std::string> branches;
  for(auto const& particle : particles.value())
    for(auto const& bf : branchFormat.value())
      branches.push_back(mysprintf(bf, particle.c_str()));

  INFO("Reading file " << inputFilename);
  INFO("Outputting file: " << outputFilename);
  TFile *f = TFile::Open(inputFilename.value().c_str(), "READ");
  INFO("Reading tree " << treeName.value());

  TTree *in_tree = (TTree *)f->Get(treeName.value().c_str());
  in_tree->SetBranchStatus("*", 1);
  for(auto &frie : friends.value() ){
      auto tokens = split(frie, ':');
      in_tree->AddFriend(tokens[1].c_str(), tokens[0].c_str());
    }

  if(inputFilename == "")
    FATAL("No input specified in options");
  if(treeName == "")
    FATAL("No tree specified in options");
  if(outputFilename == "")
    FATAL("No output specified in options");
  if(f == nullptr)
    FATAL(inputFilename.value() + " not found");
  if(in_tree == nullptr)
    FATAL(treeName.value() + " not found");

  INFO("Got tree " << inputFilename << ":" << treeName);
  std::string cut = "";
  for(auto const& c : cuts.value())
    cut += c;
  INFO("Using cut = " << cut);

  in_tree->Draw(">>elist", cut.c_str());
  TEventList *elist = (TEventList *)gDirectory->Get("elist");
  INFO("Total efficiency = " << elist->GetN() / (double)in_tree->GetEntries());

  std::vector<size_t> eventsToTake;

  if(rejectMultipleCandidates)
    {
      ULong64_t totCandidate;
      ULong64_t eventNumber;
      UInt_t runNumber;
      std::map<std::pair<ULong64_t, UInt_t>, std::vector<unsigned int>> multipleCandidateIds;

      in_tree->SetBranchStatus("*", 0);
      in_tree->SetBranchStatus("totCandidates", 1);
      in_tree->SetBranchAddress("totCandidates", &totCandidate);
      in_tree->SetBranchStatus("eventNumber", 1);
      in_tree->SetBranchAddress("eventNumber", &eventNumber);
      in_tree->SetBranchStatus("runNumber", 1);
      in_tree->SetBranchAddress("runNumber", &runNumber);

      INFO("Doing multiple candidate rejection ....");
      for(int i = 0; i < elist->GetN(); ++i)
        {
          if(i % 100000 == 0)
            INFO("Processed: " << i << " events");
          unsigned int entry = elist->GetEntry(i);
          in_tree->GetEntry(entry);
          if(totCandidate == 1)
            {
              eventsToTake.push_back(entry);
            }
          else
            {
              auto evtId = std::make_pair(eventNumber, runNumber);
              multipleCandidateIds[evtId].push_back(entry);
            }
        }

      std::mt19937 rng(7); // random-number engine used (Mersenne-Twister in this case)
      for(auto &candidate : multipleCandidateIds)
        {
          unsigned int nCand = candidate.second.size();
          unsigned int j = nCand == 1 ? 0 : std::uniform_int_distribution<int>(0, nCand - 1)(rng);
          eventsToTake.push_back(candidate.second[j]);
        }
      std::sort(eventsToTake.begin(), eventsToTake.end());
      INFO("Events before multiple candidate rejection = " << elist->GetN() << " after = " << eventsToTake.size());
    }
  else
    {
      for(int i = 0; i < elist->GetN(); ++i)
        {
          eventsToTake.push_back(elist->GetEntry(i));
        }
    }

  EventList evts(in_tree, evtType, Branches(branches), EntryList(eventsToTake), GetGenPdf(false), ApplySym(true), ExtraBranches(monitorBranches),
                 IdBranches(idBranches), InputUnits(units == "MeV" ? Units::MeV : Units::GeV));

  INFO("Branches = [" << vectorToString(branches, ", ") << "]");

  in_tree->SetBranchStatus("*", 0);
  INFO("Constructing eventList");

  if(motherID != "")
    {
      bool neg = motherID.value()[0] == '-';
      INFO("Converting " << evtType.mother() << " " << eventsToTake.size() << " " << evts.size());
      TreeReader tr(in_tree);
      int id = 0;
      tr.setBranch(neg ? motherID.value().substr(1, motherID.value().size() - 1) : motherID.value(), &id);
      for(unsigned int i = 0; i < eventsToTake.size(); ++i)
        {
          tr.getEntry(eventsToTake[i]);
          if(neg ? id > 0 : id < 0)
            invertParity(evts[i], evtType.size());
        }
    }

  INFO("Done building event list, size = " << evts.size());

  if(usePIDCalib)
    {
      INFO("Getting event weights from PID calib");
      std::string stub_path = inputFilename.value().substr(0, inputFilename.value().find_last_of('/'));
      INFO(stub_path);
      std::vector<TFile *> files;
      std::vector<TTree *> trees;
      std::vector<Float_t> weights(particles.value().size(), 0);
      for(unsigned int i = 0; i < particles.value().size(); ++i)
        {
          files.push_back(TFile::Open((stub_path + "/pidCalib_" + particles.value()[i] + "_repacked.root").c_str()));
          trees.push_back((TTree *)(*files.rbegin())->Get("CalibTool_PIDCalibTree"));
          (*trees.rbegin())->SetBranchAddress((particles.value()[i] + "_PIDCalibEffWeight").c_str(), &(weights[i]));
        }
      for(unsigned int i = 0; i < eventsToTake.size(); ++i)
        {
          double weight = 1;
          unsigned int evt = eventsToTake[i];
          for(auto &t : trees)
            {
              if(evt > t->GetEntries())
                {
                  ERROR("Accessing out of bounds data - something has gone HORRIBLY wrong");
                }
              t->GetEntry(evt);
            }
          for(auto &w : weights)
            {
              if(w < 0)
                {
                  ERROR("PID weight for event " << evt << " = " << w);
                  evts[i].print();
                }
              weight *= w;
            }
          evts[i].setWeight(weight);
        }
    }
  if(pdfLibrary != "")
    {
      INFO("Setting generator level PDF from " << pdfLibrary);
      DynamicFCN<double(const double *, const int &)> fcn(pdfLibrary, "FCN");
      for(unsigned int i = 0; i < evts.size(); ++i)
        {
          if(i % 500000 == 0)
            INFO("Set for " << i << " events");
          evts[i].setGenPdf(fcn((const real_t *)(evts[i]), 1));
        }
    }
  INFO("Writing file: " << outputFilename);
  TFile *outputFile = TFile::Open(outputFilename.value().c_str(), "RECREATE");
  INFO("Made file :-> ");

  TTree *outputTree = evts.tree("DalitzEventList");
  outputTree->Write();

  INFO("Closing file...");
  outputFile->Close();
  TFile *outputPlotFile = TFile::Open(plotsName.value().c_str(), "RECREATE");
  auto projections = evtType.defaultProjections();
  for(auto &p : projections)
    {
      p(evts)->Write();
      // p( evts, WeightFunction([](auto& evt){ return 1; }), PlotOptions::Prefix("noweight") )->Write();
    }
  for(unsigned i = 0; i != evtType.size(); ++i)
    {
      Projection p([i](auto &event) { return sqrt(event.s(i)); }, "m_" + std::to_string(i), "m_" + std::to_string(i), 100, 0, 2.5);
      p(evts)->Write();
    }
  outputPlotFile->Close();
}
