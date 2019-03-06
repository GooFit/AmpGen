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

  class RecursivePhaseSpace
  {
    private:
      struct Node 
      {
        std::string name;
        int sink = {-1};
        std::shared_ptr<RecursivePhaseSpace> decayProds = {nullptr};
        Node( const std::string& _name ) : name( _name ) {};
      };
    public:
      RecursivePhaseSpace(const EventType& type);
      RecursivePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm = gRandom);
      std::vector<Node*> getFinalStates();

      void print( const size_t& offset = 0 ) const;
      void setRandom( TRandom* rand );
      AmpGen::Event makeEvent( const size_t& cacheSize = 0 );
      size_t size() const;
      EventType eventType() const ;

    private:
      PhaseSpace   m_phsp;
      unsigned int m_totalSize;
      std::string  m_name;
      EventType    m_eventType;
      std::vector<Node> m_nodes;

  };
} // namespace AmpGen

#endif
