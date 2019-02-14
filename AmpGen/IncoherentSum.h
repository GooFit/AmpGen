#ifndef AMPGEN_FASTINCOHERENTSUM_H
#define AMPGEN_FASTINCOHERENTSUM_H

#include <complex>
#include <iomanip>
#include <string>
#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class EventType;
  class FitFraction;
  class LinearErrorPropagator;
  class MinuitParameterSet;

  class IncoherentSum : public CoherentSum
  {
    public:
      IncoherentSum( const EventType& finalStates, const AmpGen::MinuitParameterSet& mps, const std::string& prefix = "Inco" );

      double getVal( const Event& evt ) const;
      double operator()( const Event& evt ) const;
      double prob( const Event& evt ) const;
      double prob_unnormalised( const Event& evt ) const;
      void prepare();
      double norm() const;
      complex_t norm(const size_t& i, const size_t& j){ return i==j ? m_normalisations.get(i, 0) : 0; }
      complex_t norm(const size_t& i) { return m_normalisations.get(i, 0); }
      double norm( const Bilinears& norms ) const;
      std::vector<FitFraction> fitFractions( const LinearErrorPropagator& linProp );
  };
} // namespace AmpGen

#endif
