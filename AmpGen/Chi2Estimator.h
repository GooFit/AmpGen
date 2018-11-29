#ifndef AMPGEN_CHI2ESTIMATOR_H
#define AMPGEN_CHI2ESTIMATOR_H

#include <functional>
#include <string>
#include <vector>

#include "AmpGen/BinDT.h"

namespace AmpGen
{
  class EventList; 
  class EventType; 
  class Event;

  class Chi2Estimator
  {
  private: 
    double  m_chi2;
    size_t  m_nBins;
    BinDT   m_binning;

  public:
    Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
                   const std::function<double( const Event& )>& fcn, const unsigned int& minEvents = 10 );

    Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
                   const std::function<double( const Event& )>& fcn, const std::string& filename );
    double chi2() { return m_chi2; }
    double nBins() { return m_nBins; }
    void writeBinningToFile( const std::string& filename ) { m_binning.serialize( filename ); }
    void doChi2( const EventList& dataEvents, const EventList& mcEvents,
                 const std::function<double( const Event& )>& fcn );
  };
} // namespace AmpGen

#endif /* end of include guard: AMPGEN_CHI2ESTIMATOR_H */
