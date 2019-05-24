#include "AmpGen/RecursivePhaseSpace.h"

#include <cmath>

#include "AmpGen/Kinematics.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Event.h"

using namespace AmpGen;

RecursivePhaseSpace::RecursivePhaseSpace(const EventType& type) : m_eventType(type){}

RecursivePhaseSpace::RecursivePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm)
  : m_name( decayChain.name() ), m_eventType( type )
{
  std::vector<double> masses;
  m_totalSize = decayChain.getFinalStateParticles().size();
  for ( auto& d : decayChain.daughters() ) {
    masses.push_back( d->mass() );
    m_nodes.emplace_back( d->name() );
    auto decayForThisStep = std::make_shared<RecursivePhaseSpace>( *d, type );
    if ( decayForThisStep->size() != 0 ) {
      m_nodes.rbegin()->decayProds = decayForThisStep;
      INFO( "Decaying particle: " << decayChain.name() << " to " << decayChain.daughters().size() << " decay products" );
      INFO( "qsTree = " << decayChain.uniqueString() );
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
  m_phsp.setDecay( decayChain.mass(), masses );
  setRandom( rndm );
}

AmpGen::Event RecursivePhaseSpace::makeEvent( const size_t& cacheSize )
{
  AmpGen::Event evt = m_phsp.makeEvent( cacheSize );
  AmpGen::Event rt( 4 * m_eventType.size(), cacheSize );
  for (size_t i = 0; i < m_nodes.size(); ++i ) {
    auto& segment = m_nodes[i];
    double px     = evt[4*i + 0];
    double py     = evt[4*i + 1];
    double pz     = evt[4*i + 2];
    double pE     = evt[4*i + 3];
    if ( segment.decayProds == nullptr) {
      if( segment.sink != -1 ){
        rt[4*segment.sink + 0] = px;
        rt[4*segment.sink + 1] = py;
        rt[4*segment.sink + 2] = pz;
        rt[4*segment.sink + 3] = pE;
      }
    } else {
      auto evtTmp = segment.decayProds->makeEvent(cacheSize);
      double v    = sqrt( px * px + py * py + pz * pz ) / pE;
      boost( evtTmp, std::tuple<double,double,double>(px, py, pz), v );
      for(size_t j = 0; j < rt.size(); ++j) rt[j] += evtTmp[j];
    }
    evt.setGenPdf(1);
  }
  return rt;
}

std::vector<RecursivePhaseSpace::Node*> RecursivePhaseSpace::getFinalStates()
{
  std::vector<RecursivePhaseSpace::Node*> rt;
  for ( auto& f : m_nodes ) {
    if ( f.decayProds == nullptr ) rt.push_back( &f );
    else {
      auto segs = f.decayProds->getFinalStates();
      for ( auto& s : segs ) rt.push_back( s );
    }
  }
  return rt;
}

void RecursivePhaseSpace::print( const size_t& offset ) const
{
  INFO( std::string( offset, '-' ) << " Element for " << m_name << " daughter elements:" );
  for ( auto& d : m_nodes )
    if ( d.decayProds != nullptr )
      d.decayProds->print( offset + 4 );
    else
      INFO( std::string( offset + 4, '-' ) << " Final state: " << d.name << " sink = " << d.sink );
}

size_t RecursivePhaseSpace::size() const { return m_phsp.size(); }

void RecursivePhaseSpace::setRandom( TRandom* rand )
{
  m_phsp.setRandom( rand );
  for ( auto& node : m_nodes )
    if ( node.decayProds != nullptr ) node.decayProds->setRandom( rand );
}
EventType RecursivePhaseSpace::eventType() const { return m_eventType; }
