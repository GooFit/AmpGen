#include <TH1.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Factory.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/ExtendLikelihoodBase.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Property.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"

#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#if ENABLE_AVX2
#include "AmpGen/EventListSIMD.h"
using EventList_type = AmpGen::EventListSIMD;
#else
#include "AmpGen/EventList.h"
using EventList_type = AmpGen::EventList;
#endif

#include "TFile.h"
#include "TRandom3.h"

using namespace AmpGen;
using strings = std::vector<std::string>;

std::vector<ThreeBodyCalculator> threeBodyCalculators(MinuitParameterSet &mps) {
  std::vector<std::string> threeBodiesToIntegrate = Property<strings>(nullptr, "ThreeBodiesToIntegrate");
  std::vector<ThreeBodyCalculator> calculators;
  for(auto &v : threeBodiesToIntegrate) calculators.emplace_back(v, mps);
  return calculators;
}

void randomiseStartingPoint(MinuitParameterSet &MPS, TRandom3 &rand, bool SplineOnly = false) {
  double range = 5;
  for(auto &param : MPS) {
    if(!param->isFree() == 0) continue;
    if(SplineOnly && param->name().find("::Spline::") == std::string::npos) continue;
    range = param->maxInit() - param->minInit();
    param->setInit(range * rand.Rndm() + param->meanInit());
    std::cout << *param << std::endl;
  }
}

template <typename SIGPDF> void addExtendedTerms(Minimiser &mini, SIGPDF &pdf, MinuitParameterSet &mps) {
  std::vector<std::string> llConfigs = Property<std::string>(nullptr, "LLExtend");

  for(const auto &ll_config : llConfigs) {
    auto ll_name = split(ll_config, ' ')[0];
    auto ll_term = Factory<ExtendLikelihoodBase>::get(ll_name);
    if(ll_term != nullptr) {
      ll_term->configure(ll_config, pdf, mps);
      mini.addExtendedTerm(ll_term);
    } else {
      ERROR("LL term : " << ll_name << " not recognised");
    }
  }
}

template <typename PDF> FitResult *doFit(PDF &&pdf, EventList &data, EventList &mc, MinuitParameterSet &MPS) {
  INFO("Type = " << type_string<PDF>());
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time = std::clock();
  pdf.setEvents(data);

  Minimiser mini(pdf, &MPS);
  // addExtendedTerms( mini, std::get<0>( pdf.pdfs() ), MPS );
  auto threeBodyShapes = threeBodyCalculators(MPS);
  unsigned int updateWidth = Property<unsigned>(nullptr, "UpdateWidth", 0);
  unsigned int nIterations = Property<unsigned>(nullptr, "nIterations", 0);
  std::vector<std::string> SlowParams = Property<strings>(nullptr, "Release");
  bool makePlots = Property<bool>(nullptr, "MakePlots", true);

  if(updateWidth) {
    for(auto &shape : threeBodyShapes) shape.updateRunningWidth(MPS);
  }
  std::vector<MinuitProxy> slowParamPtrs;
  if(nIterations != 0) {
    for(auto &param : SlowParams) {
      auto it = MPS.find(param);
      if(it.isValid()) {
        slowParamPtrs.push_back(it);
        it->fix();
      } else {
        WARNING("Trying to release non-existent parameter: " << param);
      }
    }
  }
  INFO("Fitting PDF with " << pdf.nPDFs() << " components, iterating "
                           << " " << nIterations + 1 << " times");
  for(unsigned int iteration = 0; iteration < nIterations + 1; ++iteration) {
    mini.doFit();
    if(iteration == 0 && nIterations != 0) {
      for(auto &shape : threeBodyShapes) shape.updateRunningWidth(MPS);
      for(auto &param : slowParamPtrs) param->setFree(); /// release the parameter ///
    }
  }

  FitResult *fr = new FitResult(mini);

  if(makePlots) {
    auto ep = fr->getErrorPropagator();

    unsigned int counter = 1;
    for_each(pdf.pdfs(), [&](auto &f) {
      auto tStartIntegral2 = std::chrono::high_resolution_clock::now();
      auto mc_plot3
        = mc.makeProjections(mc.eventType().defaultProjections(100), WeightFunction(f), PlotOptions::Prefix("tMC_Category" + std::to_string(counter)));
      auto tEndIntegral2 = std::chrono::high_resolution_clock::now();
      double t2 = std::chrono::duration<double, std::milli>(tEndIntegral2 - tStartIntegral2).count();
      INFO("Time for plots = " << t2);

      for(auto &plot : mc_plot3) {
        plot->Scale((data.integral() * f.getWeight()) / plot->Integral());
        plot->Write();
      }
      counter++;
    });
  }
  Chi2Estimator chi2(data, mc, pdf, MinEvents(15));
  fr->addChi2(chi2.chi2(), chi2.nBins());

  auto twall_end = std::chrono::high_resolution_clock::now();
  double time_cpu = (std::clock() - time) / (double)CLOCKS_PER_SEC;
  double tWall = std::chrono::duration<double, std::milli>(twall_end - time_wall).count();
  INFO("Wall time = " << tWall / 1000.);
  INFO("CPU  time = " << time_cpu);
  fr->print();
  return fr;
}

int main(int argc, char *argv[]) {
  OptionsParser::setArgs(argc, argv);

  const std::string dataFile = Property<std::string>(nullptr, "DataSample", "", "Name of file containing data sample to fit.");
  const std::string mcFile = Property<std::string>(nullptr, "SimSample", "", "Name of file containing normalisation sample.");
  const std::string fmcFile = Property<std::string>(nullptr, "FlatMC", "", "Name of file containing events for computing physics integrals");
  const std::string logFile = Property<std::string>(nullptr, "LogFile", "Fitter.log", "Name of the output log file");
  const std::string plotFile = Property<std::string>(nullptr, "Plots", "plots.root", "Name of the output plot file");

  [[maybe_unused]] const size_t nThreads = Property<size_t>(nullptr, "nCores", 8, "Number of threads to use");

  const size_t NBins = Property<size_t>(nullptr, "nBins", 100, "Number of bins used for plotting.");
  const bool perturb = Property<bool>(nullptr, "Perturb", 0, "Flag to randomise starting parameters.");
  const size_t seed = Property<size_t>(nullptr, "Seed", 0, "Random seed used");

  std::vector<std::string> evtType_particles
    = Property<strings>(nullptr, "EventType", {}, "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m");

  const std::string cut = Property<std::string>(nullptr, "Cut", "1");
  const std::string simCut = Property<std::string>(nullptr, "SimCut", "1");
  bool BAR = Property<bool>(nullptr, "Bar", false);
  const std::string units = Property<std::string>(nullptr, "Units", "GeV");

  INFO("Output : " << logFile << " plots = " << plotFile);

  TRandom3 rndm;
  rndm.SetSeed(seed);
  gRandom = &rndm;

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
  INFO("Setting " << nThreads << " fixed threads for OpenMP");
  omp_set_dynamic(0);
#endif

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if(dataFile == "") {
    ERROR("No input data selected");
    return -1;
  }
  if(mcFile == "") WARNING("No input simulation selected; using PHSP sample");

  if(perturb) {
    for(auto &param : MPS) {
      if(!param->isFree()) continue;
      param->setCurrentFitVal(rndm.Gaus(param->mean(), param->err()));
    }
  }
  EventType evtType(evtType_particles);
  INFO("Signal  = " << evtType << " OS = " << evtType.conj(true));
  CoherentSum sig(evtType, MPS);
  IncoherentSum bkg(evtType, MPS, "Inco");
  CoherentSum misID(evtType.conj(true), MPS);

  INFO("fPDF = " << MPS["fPDF"]->mean() << " "
                 << "fComb = " << MPS["fComb"]->mean() << " "
                 << "fMisID = " << MPS["fMisID"]->mean());

  sig.setWeight(MPS["fPDF"]);
  bkg.setWeight(MPS["fComb"]);
  misID.setWeight(MPS["fMisID"]);

  EventList events(dataFile, !BAR ? evtType : evtType.conj(), Filter(cut));
  EventList eventsMC = mcFile == "" ? EventList(evtType) : EventList(mcFile, !BAR ? evtType : evtType.conj(), Filter(simCut));

  auto scale_transform = [](auto &event) {
    for(size_t x = 0; x < event.size(); ++x) event[x] /= 1000.;
  };
  if(units == "Mev") {
    INFO("Changing units from MeV -> GeV");
    events.transform(scale_transform);
  }
  eventsMC.transform(scale_transform);

  INFO("Data events: " << events.size());
  INFO("MC events  : " << eventsMC.size());
  if(mcFile == "") {
    eventsMC = Generator<>(evtType, &rndm).generate(5e6);
    INFO("Generated: " << eventsMC.size() << " events for integrals");
  }

  sig.setMC(eventsMC);
  bkg.setMC(eventsMC);
  misID.setMC(eventsMC);

  TFile *output = TFile::Open(plotFile.c_str(), "RECREATE");
  output->cd();
  FitResult *fr = nullptr;

  if(MPS["fPDF"]->mean() == 1.0)
    fr = doFit(make_pdf(sig), events, eventsMC, MPS);

  else if(MPS["fPDF"]->mean() == 0.0 && MPS["fComb"]->mean() == 1.0)
    fr = doFit(make_pdf(bkg), events, eventsMC, MPS);

  else if(MPS["fMisID"]->mean() == 0.0)
    fr = doFit(make_pdf(sig, bkg), events, eventsMC, MPS);

  else if(MPS["fMisID"]->mean() != 0.0)
    fr = doFit(make_pdf(sig, bkg, misID), events, eventsMC, MPS);

  if(fr == nullptr) {
    ERROR("Fit fails");
    return -1;
  }

  INFO("Completed fit; calculating additional observables");

  //  EventList* fmc               = &eventsMC;
  //  if ( fmcFile != mcFile ) fmc = new EventList( fmcFile, evtType );
  //  EventList& flatMC            = *fmc;
  //  sig.reset( true ); //// reset PDFs to ensure correct cache state
  //  sig.setMC( flatMC );
  //  sig.prepare();
  if(MPS["fComb"]->mean() != 1) fr->addFractions(sig.fitFractions(fr->getErrorPropagator()));
  //  else
  //    fr->addFractions( bkg.fitFractions( fr->getErrorPropagator() ) );

  fr->writeToFile(logFile);
  output->cd();
  auto plots = events.makeDefaultProjections(PlotOptions::Prefix("Data_"), PlotOptions::Bins(NBins));
  for(auto &plot : plots) plot->Write();

  output->Write();
  output->Close();

  INFO("Finalising output");
  return 0;
}
