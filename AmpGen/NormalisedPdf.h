#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/MinuitParameterSet.h"
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
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ErrorPropagator.h"

#include <TFile.h>
#include <TGraph2D.h>
#include <TRandom3.h>
#include <TRandom.h>
#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#include "AmpGen/DTEvent.h"
#include "AmpGen/FixedLibPdf.h"
#include "AmpGen/DTYieldCalculator.h"

using namespace AmpGen;
using namespace std::complex_literals;
template <class PDF> struct normalised_pdf {
  PDF       pdf; 
  complex_t norm;
  normalised_pdf() = default; 
  normalised_pdf(const EventType& type, MinuitParameterSet& mps, const DTYieldCalculator& yc) : pdf(type, mps)
  {
    ProfileClock pc; 
    auto normEvents = Generator<PhaseSpace>(type).generate(1e6);
    double n =0;
    pdf.prepare();
    #pragma omp parallel for reduction(+:n)
    for(size_t i =0; i<normEvents.size(); ++i) 
      n += std::norm(pdf.getValNoCache(normEvents[i]));
    auto it = mps.find( type.decayDescriptor() + "::strongPhase");
    norm = sqrt(yc.bf(type)/n);
    if( it != nullptr ) norm *= exp( 1i * it->mean() * M_PI/180. );
    pc.stop();
    INFO(type << " Time to construct: " << pc << "[ms], norm = " << norm  << " " << typeof<PDF>() );
  }
  complex_t operator()(const Event& event){ return norm * pdf.getValNoCache(event); }
};

struct ModelStore {
  MinuitParameterSet*                                 mps;
  DTYieldCalculator                                   yieldCalculator;
  std::map<std::string, normalised_pdf<CoherentSum>>  genericModels;
  std::map<std::string, normalised_pdf<FixedLibPdf>>  flibModels;
  ModelStore(MinuitParameterSet* mps, const DTYieldCalculator& yc) : mps(mps), yieldCalculator(yc) {}
  template <class T> normalised_pdf<T>& get(const EventType& type, std::map<std::string, normalised_pdf<T>>& container)
  {
    auto key = type.decayDescriptor();
    if( container.count(key) == 0 ) 
      container[key] = normalised_pdf<T>(type, *mps, yieldCalculator);
    return container[key]; 
  }
  template <class T>    normalised_pdf<T>& find(const EventType& type);
};