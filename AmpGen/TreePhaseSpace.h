#ifndef AMPGEN_TREEPHASESPACE_H
#define AMPGEN_TREEPHASESPACE_H

#include <memory.h>
#include <stddef.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include <random>

#include "AmpGen/EventType.h"
#include "AmpGen/Particle.h"
#include "AmpGen/DiscreteDistribution.h"
#include <TRandom3.h>

namespace AmpGen
{
  class Particle;
  class Event; 
  class DecayChainStackBase; 
  class EventListSIMD; 
  /** @class TreePhaseSpace
    @brief Generator of events where the phase space is decomposed into a series of subtrees.  
    @decription Generates events using the decomposition of the phase space 
    into a series of two-body phase spaces, where the invariant mass of the two-body space 
    can be generated according to, for example, a simplified Breit-Wigner model. 
    This can then be used to generate a more complex model using accept-reject, which 
    allows for a much more efficient generation of narrow peaks. 

    Currently, each of the available channels is given the same weight relative to phase space, 
    ideally feedback could be given from the generator phase to focus on the more efficient channels, 
    i.e. those that have larger contributions to the full amplitude.   
    */
  class TreePhaseSpace
  {
    public: 
      explicit TreePhaseSpace(const EventType& type);
      TreePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm = nullptr );
      TreePhaseSpace(const std::vector<Particle>& decayChains, const EventType& type, TRandom* rndm = nullptr);
      ~TreePhaseSpace(); 

      void setRandom( TRandom* rand );
      Event makeEvent();
      size_t size() const;
      EventType eventType() const ;
      const DecayChainStackBase* operator[](const unsigned i) const { return m_gen[i]; }      
      void recalculate_weights( const EventListSIMD& events); 
    private:
      void initialise_weights(); 
      std::vector<DecayChainStackBase*> m_gen; 
      TRandom3*             m_rand   {nullptr}; 
      EventType             m_type;         ///< EventType to generate
      DiscreteDistribution  m_dice;  
      std::vector<unsigned> m_generatorRecord;  
      double                m_wmax   {0}; 
  };
} // namespace AmpGen

#endif
