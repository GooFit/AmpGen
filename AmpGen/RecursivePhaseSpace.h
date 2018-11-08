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

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/PhaseSpace.h"
#include "TRandom3.h"
#include "AmpGen/Event.h"

class TRandom;

namespace AmpGen
{
  class RecursivePhaseSpace
  {
  private:
    struct RecursivePhaseSpaceNode {
      std::string name;
      int sink;
      std::shared_ptr<RecursivePhaseSpace> decayProds;
      RecursivePhaseSpaceNode( const std::string& _name ) : name( _name ), sink( -1 ), decayProds( nullptr ) {}
    };
    PhaseSpace m_phsp;
    unsigned int m_totalSize;
    std::string m_name;
    EventType m_eventType;
    std::vector<RecursivePhaseSpaceNode> m_nodes;

  public:
    RecursivePhaseSpace( const EventType& type ) : m_eventType( type ) {}

    void SetDecay( const double& mother, const std::vector<double>& daughters ) {}

    RecursivePhaseSpace( const Particle& decayChain, const EventType& type, TRandom*  rndm = gRandom )
        : m_name( decayChain.name() ), m_eventType( type )
    {
      std::vector<double> masses;
      m_totalSize = decayChain.getFinalStateParticles().size();
      DEBUG( "Decaying particle: " << decayChain.name() << " to " << decayChain.daughters().size()
                                   << " decay products" );
      DEBUG( "qsTree = " << decayChain.uniqueString() );
      for ( auto& d : decayChain.daughters() ) {
        masses.push_back( d->mass() );
        m_nodes.emplace_back( d->name() );
        auto decayForThisStep = std::make_shared<RecursivePhaseSpace>( *d, type );

        if ( decayForThisStep->size() != 0 ) {
          m_nodes.rbegin()->decayProds = decayForThisStep;
        }
      }
      auto fs = getFinalStates();
      std::vector<int> used( m_eventType.size(), 0 );
      for ( auto& f : fs ) {
        for ( unsigned int i = 0; i < m_eventType.size(); ++i ) {
          if ( !used[i] && f->name == m_eventType[i] ) {
            f->sink = i;
            used[i] = true;
            break;
          }
        }
      }
      m_phsp.SetDecay( decayChain.mass(), masses );
      setRandom( rndm );
    }
    std::vector<RecursivePhaseSpaceNode*> getFinalStates()
    {
      std::vector<RecursivePhaseSpaceNode*> rt;
      for ( auto& f : m_nodes ) {
        if ( f.decayProds == nullptr )
          rt.push_back( &f );
        else {
          auto segs = f.decayProds->getFinalStates();
          for ( auto& s : segs ) rt.push_back( s );
        }
      }
      return rt;

    }

    void print( const size_t offset = 0 ) const
    {
      INFO( std::string( offset, '-' ) << " Element for " << m_name << " daughter elements:" );
      for ( auto& d : m_nodes )
        if ( d.decayProds != nullptr )
          d.decayProds->print( offset + 4 );
        else
          INFO( std::string( offset + 4, '-' ) << " Final state: " << d.name << " sink = " << d.sink );
    }
    AmpGen::Event makeEvent( const size_t& cacheSize = 0 );
    size_t size() const { return m_phsp.size(); }
    void setRandom( TRandom* rand )
    {
      m_phsp.setRandom( rand );
      for ( auto& node : m_nodes )
        if ( node.decayProds != nullptr ) node.decayProds->setRandom( rand );
    }
    EventType eventType() const { return m_eventType; }
  };
} // namespace AmpGen

#endif
