#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>

#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Property.h"
#include "AmpGen/PolarisedSum.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AmpGen;
using strings = std::vector<std::string>;

void invertParity(Event &event, const size_t &nParticles) {
  for(size_t i = 0; i < nParticles; ++i) {
    event[4 * i + 0] = -event[4 * i + 0];
    event[4 * i + 1] = -event[4 * i + 1];
    event[4 * i + 2] = -event[4 * i + 2];
  }
}

void randomBoost(Event &event, TRandom3 *rndm) {
  auto v = std::make_tuple(rndm->Uniform(), rndm->Uniform(), rndm->Uniform());
  auto beta = rndm->Uniform();
  Event pol(4);
  pol[0] = -0.104734;
  pol[1] = -0.328773;
  pol[2] = 0.617203;
  pol[3] = 0.0;
  boost(pol, v, beta);
  pol.print();
  boost(event, v, beta);
}

void randomRotation(Event &event, TRandom3 *rndm) { rotate(event, std::make_tuple(rndm->Uniform(), rndm->Uniform(), rndm->Uniform()), rndm->Uniform()); }

template <typename fcn_type> void writeRefFile(const std::string &filename, fcn_type &fcn, const EventList &events) {
  using json = nlohmann::json;
  json output;
  for(unsigned i = 0; i != events.size(); ++i) {
    output["event_" + std::to_string(i)] = events[i].data();
    output["pdf_" + std::to_string(i)] = fcn(events[i]);
    for(const auto &elem : fcn.matrixElements()) {
      auto indices = fcn.cache().find(elem.name());
      std::vector<double> this_cache;
      for(const auto &index : indices) {
        auto v = utils::at(fcn.cache()(i / utils::size<real_v>::value, index), i % utils::size<real_v>::value);
        this_cache.emplace_back(std::real(v));
        this_cache.emplace_back(std::imag(v));
      }
      output[elem.name() + "_" + std::to_string(i)] = this_cache;
    }
  }
  std::ofstream os(filename);
  os << output.dump(4) << std::endl;
  os.close();
}

template <typename FCN> void debug(FCN &sig, EventList &accepted) {
  INFO("Debugging: ");
  unsigned eventToDebug = 0;
  sig.setEvents(accepted);
  sig.prepare();
  accepted[eventToDebug].print();
  sig.debug(accepted[eventToDebug]);

  /*
  INFO("Parity: " );

  for( unsigned int i = 0 ; i != accepted.size(); ++i )
    invertParity(accepted[i], accepted.eventType().size() );
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  */
  INFO("Random rotation:");
  auto old_event = accepted[eventToDebug];
  randomRotation(accepted[eventToDebug], new TRandom3(eventToDebug));
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug(accepted[eventToDebug]);
  accepted[eventToDebug] = old_event;
  INFO("Random boost:");

  randomBoost(accepted[eventToDebug], new TRandom3(eventToDebug));
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug(accepted[eventToDebug]);
}

int main(int argc, char **argv) {
  OptionsParser::setArgs(argc, argv);

  EventType eventType(Property<strings>(nullptr, "EventType", {}, "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m"),
                      Property<bool>(nullptr, "GenerateTimeDependent", false, "Flag to include possible time dependence of the amplitude"));

  int seed = Property<int>(nullptr, "Seed", 156);
  std::string infile = Property<std::string>(nullptr, "InputFile", "");
  std::string refFileOutput = Property<std::string>(nullptr, "RefFileOutput", "");
  std::string input_units = Property<std::string>(nullptr, "Units", "GeV");
  std::string type = Property<std::string>(nullptr, "Type", "CoherentSum");

  std::vector<double> event = Property<std::vector<double>>(nullptr, "Event", {});
  bool conj = Property<bool>(nullptr, "conj", false);
  bool add_conj = Property<bool>(nullptr, "AddConj", false);

  bool verbose = Property<bool>(nullptr, "CoherentSum::Debug", false) || Property<bool>(nullptr, "PolarisedSum::Debug", false);

  INFO("Using verbose mode: " << verbose);
  AmpGen::MinuitParameterSet MPS;
  MPS.loadFromStream();

  TRandom3 *rndm = new TRandom3(seed);

  if(conj) {
    eventType = eventType.conj();
    INFO(eventType);
    AddCPConjugate(MPS);
  } else if(add_conj) {
    AddCPConjugate(MPS);
  }
  INFO("EventType = " << eventType);

  EventList accepted = infile == "" ? EventList(eventType) : EventList(infile, eventType);

  if(input_units == "MeV" && infile != "")
    accepted.transform([](auto &event) {
      for(unsigned i = 0; i < event.size(); ++i) event[i] /= 1000;
    });
  if(infile == "") {
    accepted = Generator<PhaseSpace>(eventType, rndm).generate(16);
    for(unsigned i = 0; i != 16; ++i) accepted[i].setIndex(i);
  }
  if(event.size() != 0) accepted[0].set(event.data());

  if(type == "PolarisedSum") {
    PolarisedSum sig(eventType, MPS);
    sig.setEvents(accepted);
    sig.prepare();
    if(refFileOutput != "") writeRefFile(refFileOutput, sig, accepted);
    debug(sig, accepted);
    sig.setMC(accepted);
    INFO("norm = " << sig.norm());
  } else if(type == "CoherentSum") {
    CoherentSum sig(eventType, MPS);
    sig.setEvents(accepted);
    sig.prepare();
    if(refFileOutput != "") writeRefFile(refFileOutput, sig, accepted);
    debug(sig, accepted);
    INFO("A(x) = " << sig.getValNoCache(accepted[0]));
  } else if(type == "IncoherentSum") {
    IncoherentSum sig(eventType, MPS, "Inco");
    sig.setMC(accepted);
    sig.prepare();
    debug(sig, accepted);
    INFO("norm = " << sig.norm());
  } else {
    ERROR("Type: " << type << " is not recognised");
  }
}
