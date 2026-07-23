#include <Rtypes.h>
#include <TH1.h>
#include <dlfcn.h>
#include <memory>
#include <string>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

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
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/enum.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/Property.h"
#include "AmpGen/ApplicationOptions.h"

#if ENABLE_AVX
#include "AmpGen/EventListSIMD.h"
using EventList_t = AmpGen::EventListSIMD;
#else
#include "AmpGen/EventList.h"
using EventList_t = AmpGen::EventList;
#endif

using namespace AmpGen;

namespace AmpGen {
  make_enum(pdfTypes, CoherentSum, IncoherentSum, PolarisedSum) make_enum(phspTypes, PhaseSpace, RecursivePhaseSpace, TreePhaseSpace)
}

template <class T> void generateSource(T &pdf, const std::string &sourceFile, MinuitParameterSet &mps) {
  bool normalise = Property<bool>(nullptr, "Normalise", true);
  double safetyFactor = Property<double>(nullptr, "SafetyFactor", 3);
  int seed = Property<int>(nullptr, "Seed", 1);
  size_t nEvents = Property<size_t>(nullptr, "NormEvents", 1000000);

  double norm = 1;
  if(normalise) {
    INFO("Normalising PDF");
    TRandom3 rnd(seed);
    unsigned d_i = pdf.eventType().dim().first;
    Generator<PhaseSpace> phsp(pdf.eventType());
    phsp.setRandom(&rnd);
    EventList_t normEvents = phsp.generate(nEvents);
    if constexpr(std::is_same<T, CoherentSum>::value) pdf.prepare();
    double pMax = 0;
    for(auto &evt : normEvents) {
      if constexpr(std::is_same<T, PolarisedSum>::value) {
        if(d_i > 1) {
          double px, py, pz;
          rnd.Sphere(px, py, pz, rnd.Uniform(0, 1));
          mps["Px"]->setCurrentFitVal(px);
          mps["Py"]->setCurrentFitVal(py);
          mps["Pz"]->setCurrentFitVal(pz);
          pdf.transferParameters();
        }
      }
      double n = 0;
      if constexpr(std::is_same<T, CoherentSum>::value) n = std::norm(pdf.getValNoCache(evt));
      if constexpr(std::is_same<T, PolarisedSum>::value) n = pdf.getValNoCache(evt);
      if(n > pMax) pMax = n;
    }
    norm = pMax * safetyFactor;
    INFO("Making binary with " << pMax << " x safety factor = " << safetyFactor);
  }
  mps.resetToInit();
  pdf.generateSourceCode(sourceFile, norm, true);
}

template <typename pdf_t> Particle getTopology(const pdf_t &pdf) { return pdf.matrixElements()[0].decayTree.quasiStableTree(); }

template <typename pdf_t> std::vector<Particle> getDecayChains(const pdf_t &pdf) {
  std::vector<Particle> channels;
  for(auto &chain : pdf.matrixElements()) channels.push_back(chain.decayTree);
  return channels;
}

template <typename pdf_t>
void generateEvents(EventList &events, pdf_t &pdf, const phspTypes &phsp_type, const size_t &nEvents, const size_t &blockSize, TRandom *rndm,
                    const bool &normalise = true) {
  auto fill = [&](auto &generator) mutable {
    generator.setRandom(rndm);
    generator.setBlockSize(blockSize);
    generator.setNormFlag(normalise);
    generator.fillEventList(pdf, events, nEvents);
  };

  if(phsp_type == phspTypes::PhaseSpace) {
    Generator<PhaseSpace> signalGenerator(events.eventType());
    fill(signalGenerator);
  } else if(phsp_type == phspTypes::RecursivePhaseSpace) {
    Generator<RecursivePhaseSpace> signalGenerator(getTopology(pdf), events.eventType());
    fill(signalGenerator);
  } else if(phsp_type == phspTypes::TreePhaseSpace) {
    Generator<TreePhaseSpace> signalGenerator(getDecayChains(pdf), events.eventType());
    fill(signalGenerator);
  } else {
    FATAL("Phase space configuration: " << phsp_type << " is not supported");
  }
}

int main(int argc, char **argv) {
  OptionsParser::setArgs(argc, argv);
  ApplicationOptions options;

  std::string phspType_hs = helpStringOptions(
    "Phase-space generator to use:", std::make_pair(phspTypes::PhaseSpace, "Phase space generation based on Raubold-Lynch algorithm (recommended).\0"),
    std::make_pair(phspTypes::TreePhaseSpace,
                   "Divides the phase-space into a series of quasi two-body phase-spaces for efficiently generating narrow states.\0"),
    std::make_pair(phspTypes::RecursivePhaseSpace,
                   "Includes possible quasi-stable particles and the phase spaces of their decay products, such as Λ baryons.\0"));

  std::string pdfType_hs
    = helpStringOptions("Type of PDF to use:", std::make_pair(pdfTypes::CoherentSum, "Describes decays of a (pseudo)scalar particle to N pseudoscalars"),
                        std::make_pair(pdfTypes::IncoherentSum, "Describes background-like contribution to pseudoscalar decay processes."),
                        std::make_pair(pdfTypes::PolarisedSum, "Describes the decay of a particle with spin to N particles carrying spin."));

  Property<std::string> decay{nullptr, "Decay", "", "Single decay written on the command line, overwrites all other options."};
  Property<size_t> nEvents{nullptr, "nEvents", 1, "Total number of events to generate"};
  Property<size_t> blockSize{nullptr, "BlockSize", 5000000, "Number of events to generate per block"};
  Property<pdfTypes> pdfType{nullptr, "Type", pdfTypes::CoherentSum, pdfType_hs};
  Property<phspTypes> phspType{nullptr, "PhaseSpace", phspTypes::PhaseSpace, phspType_hs};
  auto ext = *split(options.outfile, '.').rbegin();
  bool sourceOnly = ext == "root" ? false : true;
  Property<bool> poissonVaryYield{nullptr, "PoissonYield", false, "Vary the number of events generated by a poisson distribution"};
  Property<bool> generateTD{nullptr, "GenerateTimeDependent", false, "Flag to include possible time dependence of the amplitude"};

  Property<size_t> nBins{nullptr, "nBins", 100, "Number of bins for monitoring plots."};
  Property<std::string> observable{nullptr, "Observable", "mass2", "Observable to use for making plots {mass, mass2}"};
  Property<bool> make2DPlots{nullptr, "Make2Dplots", false, "Make two-dimensional projections."};

#ifdef _OPENMP
  omp_set_num_threads(options.nCores);
  omp_set_dynamic(0);
#endif

  TRandom3 rand;
  rand.SetSeed(options.seed + 934534);

  if(poissonVaryYield) nEvents.set(rand.Poisson(nEvents));
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  if(OptionsParser::printHelp()) return 0;

  EventType eventType;
  if(decay != "") {
    Particle p(decay);
    eventType = p.eventType();
    MPS.add(p.decayDescriptor() + "_Re", Flag::Fix, 1., 0);
    MPS.add(p.decayDescriptor() + "_Im", Flag::Fix, 0., 0);
  } else
    eventType = EventType(options.eventType_s, generateTD);

  if(options.conj) eventType = eventType.conj();
  if(options.conj || options.addCPConjugate) AddCPConjugate(MPS, options.forbidCP);

  INFO("Writing output: " << options.outfile);

  auto [dim_i, dim_f] = eventType.dim();
  if((dim_i != 1 || dim_f != 1) && pdfType == pdfTypes::CoherentSum) {
    WARNING("Either the initial or final state involves a particle that carries spin, switching to use PolarisedSum");
    pdfType.set(pdfTypes::PolarisedSum);
  }

  if(sourceOnly) {
    if(pdfType == pdfTypes::CoherentSum) {
      CoherentSum pdf(eventType, MPS);
      generateSource(pdf, options.outfile, MPS);
    }
    if(pdfType == pdfTypes::IncoherentSum) {
      CoherentSum pdf(eventType, MPS, "Inco");
      generateSource(pdf, options.outfile, MPS);
    } else if(pdfType == pdfTypes::PolarisedSum) {
      PolarisedSum pdf(eventType, MPS);
      generateSource(pdf, options.outfile, MPS);
    }
    return 0;
  }
  INFO("Generating time-dependence? " << eventType.isTimeDependent());
  EventList accepted(eventType);

  INFO("Generating events with type = " << eventType);

  if(pdfType == pdfTypes::CoherentSum) {
    CoherentSum pdf(eventType, MPS);
    if(pdf.size() == 0) FATAL("Requested model has no amplitudes");
    generateEvents(accepted, pdf, phspType, nEvents, blockSize, &rand);
  } else if(pdfType == pdfTypes::IncoherentSum) {
    CoherentSum pdf(eventType, MPS, "Inco");
    if(pdf.size() == 0) FATAL("Requested model has no amplitudes");
    generateEvents(accepted, pdf, phspType, nEvents, blockSize, &rand);
  } else if(pdfType == pdfTypes::PolarisedSum) {
    PolarisedSum pdf(eventType, MPS);
    if(pdf.size() == 0) FATAL("Requested model has no amplitudes");
    generateEvents(accepted, pdf, phspType, nEvents, blockSize, &rand);
  } else {
    FATAL("Did not recognise configuration: " << pdfType);
  }
  if(accepted.size() == 0) return -1;
  TFile *f = TFile::Open(options.outfile.value().c_str(), "RECREATE");
  accepted.tree("DalitzEventList")->Write();
  auto proj = eventType.defaultProjections(nBins, observable);
  auto plots = accepted.makeProjections(proj, PlotOptions::LineColor(kBlack));
  for(auto &plot : plots) plot->Write();
  if(make2DPlots) {
    for(size_t i = 0; i < proj.size(); ++i) {
      for(size_t j = i + 1; j < proj.size(); ++j) { accepted.makeProjection(Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack))->Write(); }
    }
  }
  INFO("Writing output file ");

  f->Close();
}
