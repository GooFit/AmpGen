#ifndef AMPGEN_INCOHERENTSUM_H
#define AMPGEN_INCOHERENTSUM_H

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
  
  /** @class IncoherentSum
      @brief An incoherent sum of resonances, traditionally considered as an approximate background model. 
      An incoherent sum of resonant contributions, where at a positiion in phase-space @f$\psi@f$, 
      the (unnormalised) probability density  is given by
      @f[ 
        \mathcal{P}(\psi) = \sum_{i} \left| g_i \mathcal{A}_i(\psi) \right|^2,    
      @f]
      where @f$\mathcal{P}(\psi)@f$ is the probability, @f$g_i@f$ is the coupling to an isobar channel, 
      and @f$\mathcal{A}_i(\psi)@f$ is the amplitude of the ith channel. 
      I dont know whats going on?
      */ 
  class IncoherentSum : public CoherentSum
  {
    public:
      /** Constructs an incoherentSum from the type of event this is expected to describe, 
         and a set of parameters that are expected to be able to describe it.  
          @param eventType The type of event that this PDF should describe
          @param mps       The parameter set of couplings, masses, etc.
          @param prefix    Prefix required for the ``head'' decays for this event type. */ 
      IncoherentSum(const EventType& eventType, const AmpGen::MinuitParameterSet& mps, const std::string& prefix = "Inco");
      
      /// Evaluates the normalised probability for an event.
      real_t prob( const Event& evt ) const;
      real_t operator()(const Event& evt) const { return prob(evt) ; }

      /// Calculates the unnormalised probability for an event. 
      real_t prob_unnormalised( const Event& evt ) const;

      /** Returns the normalisation for this PDF, given by
          @f[
            \mathcal{N} = \int d\psi \varepsilon(\psi) \mathcal{P}(\psi) \approx \sum_i \frac{\mathcal{P}(\psi_i)}{\mathcal{P}^\prime(\psi_i)} 
          @f]
          where the sum is over a simulated sample, generated with PDF @f$\mathcal{P}^\prime(\psi)@f$. */
      real_t norm() const;
      complex_t norm(const size_t& i, const size_t& j){ return i==j ? m_normalisations.get(i, 0) : 0; }
      complex_t norm(const size_t& i) { return m_normalisations.get(i, 0); }
      real_t norm( const Bilinears& norms ) const;
      std::vector<FitFraction> fitFractions( const LinearErrorPropagator& linProp );
      
      void prepare();
  };
} // namespace AmpGen

#endif
