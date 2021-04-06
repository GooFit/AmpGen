#ifndef AMPGEN_CHI2ESTIMATOR_H
#define AMPGEN_CHI2ESTIMATOR_H

#include <functional>
#include <string>
#include <vector>

#include "AmpGen/BinDT.h"

#if ENABLE_AVX
#include "AmpGen/EventListSIMD.h"
#else 
#include "AmpGen/EventList.h"
#endif


namespace AmpGen
{
  class EventType; 
  class Event;

  class Chi2Estimator
  {
    #if ENABLE_AVX
    typedef EventListSIMD EventList_type;
    #else
    typedef EventList     EventList_type;
    #endif
    public:
    template <typename... argument_types> 
    Chi2Estimator( const EventList_type& dataEvents, const EventList_type& mcEvents,
                   const std::function<double( const Event& )>& fcn, 
                   const argument_types&... args ) : m_binning(dataEvents.begin(), dataEvents.end(), ArgumentPack(args...) ) 
    {
      doChi2(dataEvents, mcEvents, fcn);
    }

    double chi2() const;
    double nBins() const; 
    void writeBinningToFile( const std::string& filename ); 
    void doChi2( const EventList_type& dataEvents, const EventList_type& mcEvents,
                 const std::function<double( const Event& )>& fcn );
  private: 
    double  m_chi2;
    size_t  m_nBins;
    BinDT   m_binning;
  };
} // namespace AmpGen

#endif /* end of include guard: AMPGEN_CHI2ESTIMATOR_H */
