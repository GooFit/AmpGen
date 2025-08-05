#include "AmpGen/Configurable.h"
#include "AmpGen/Property.h"

namespace AmpGen {

  /// default
  struct ApplicationOptions : public Configurable<ApplicationOptions> {
    virtual ~ApplicationOptions() = default;

    using strings = std::vector<std::string>;

    Property<strings> eventType_s{this, "EventType", {}, "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m"};
    Property<int> seed{this, "Seed", 0, "Random seed used in event Generation. Should always be set for pseudoexperiment generation."};
    Property<std::string> outfile{this, "Output", "Output.root", "Name of output file"};
    Property<bool> conj{this, "Conj", false, "Flag to generate the CP conjugate amplitude under the assumption of CP conservation"};
    Property<bool> addCPConjugate{this, "AddConj", false, "Flag to add all of the CP conjugate amplitudes, under the assumption of CP conservation"};
    Property<strings> forbidCP{this, "ForbidCP", {}, "Parameters to forbid the addition of corresponding CP conjugate"};
#ifdef _OPENMP
    Property<unsigned> nCores{this, "nCores", std::thread::hardware_concurrency(), "Number of cores to use (OpenMP only)"};
#endif
  };

};
