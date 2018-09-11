#include "AmpGen/RecursivePhaseSpace.h"

#include <cmath>

#include "AmpGen/Kinematics.h"

using namespace AmpGen;

AmpGen::Event RecursivePhaseSpace::makeEvent( const size_t& cacheSize )
{

  DEBUG( "Making event for decay of " << m_name << " to " << m_phsp.size()
                                      << " daughters [type size = " << m_eventType.size() << "]" );
  AmpGen::Event evt = m_phsp.makeEvent( cacheSize );
  AmpGen::Event rt( 4 * m_eventType.size(), cacheSize );
  for ( unsigned int i = 0; i < m_nodes.size(); ++i ) {
    auto& segment = m_nodes[i];
    double px     = evt[4 * i + 0];
    double py     = evt[4 * i + 1];
    double pz     = evt[4 * i + 2];
    double pE     = evt[4 * i + 3];
    if ( segment.decayProds == nullptr ) {
      rt[4 * segment.sink + 0] = px;
      rt[4 * segment.sink + 1] = py;
      rt[4 * segment.sink + 2] = pz;
      rt[4 * segment.sink + 3] = pE;
    } else {
      auto evtTmp = segment.decayProds->makeEvent( cacheSize );
      double v    = -sqrt( px * px + py * py + pz * pz ) / pE;
      boost( evtTmp, {px, py, pz}, v );
      for ( unsigned int i = 0; i < rt.size(); ++i ) rt[i] += evtTmp[i];
    }
    evt.setGenPdf( 1 );
  }
  return rt;
}
