#ifndef AMPGEN_RECURSIVEPHASESPACE_H
#define AMPGEN_RECURSIVEPHASESPACE_H

#include <memory.h>
#include <stddef.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"

#include <TRandom3.h>

class TRandom;

namespace AmpGen
{
  class Particle;
  class Event; 
  /** @class RecursivePhaseSpace
      @brief Generator of events with more complex topologies.
      @decription Generates events in the phase space of a decay with one-or-more quasistable products 
      that are further decayed, where these quasistable products are suitably long-lived 
      that their propagators can be assumed to be delta-functions, i.e. particles such as 
      @f$ \Lambda^0 , {K_\rm{S}}^0 @f$ etc. Such a process has a natural factorisation into 
      each of these steps. For each of the decays of a quasistable particle, the PhaseSpace 
      generator is used, which is itself based heavily on TGenPhaseSpace.  
  */
  class RecursivePhaseSpace
  {
    private:
      struct Node 
      {
        std::string name = {""};
        int sink         = {-1};
        std::shared_ptr<RecursivePhaseSpace> decayProds = {nullptr};
        explicit Node( const std::string& _name ) : name( _name ) {};
      };
    public:
      explicit RecursivePhaseSpace(const EventType& type);
      RecursivePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm = gRandom);
      std::vector<Node*> getFinalStates();

      void print( const size_t& offset = 0 ) const;
      void setRandom( TRandom* rand );
      Event makeEvent();
      size_t size() const;
      EventType eventType() const ;

      void provideEfficiencyReport(const std::vector<bool>& report){}
    private:
      PhaseSpace   m_phsp;
      unsigned     m_totalSize = {0};
      std::string  m_name      = {""};
      EventType    m_eventType;
      std::vector<Node> m_nodes;
  };
} // namespace AmpGen

#endif
