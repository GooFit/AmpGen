#ifndef AMPGEN_CHI2ESTIMATOR_H
#define AMPGEN_CHI2ESTIMATOR_H

#include <functional>
#include <string>
#include <vector>

#include "AmpGen/BinDT.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"

namespace AmpGen
{
  class Chi2Estimator
  {

    double m_chi2;
    unsigned int m_nBins;
    BinDT m_binning;

  public:
    double chi2() { return m_chi2; }
    double nBins() { return m_nBins; }
    struct Moment {
      double x;
      double xx;
      double N;
      std::vector<double> values;
      Moment() : x( 0 ), xx( 0 ), N( 0 ) {}
      void add( const double& value )
      {
        x += value;
        xx += value * value;
        N++;
        values.push_back( value );
      }
      void rescale( const double& val )
      {
        x *= val;
        xx *= ( val * val );
      }
      double val() { return x; }
      double var() { return N == 0 ? 0 : xx; }
    };

    void writeBinningToFile( const std::string& filename ) { m_binning.serialize( filename ); }
    void doChi2( const EventList& dataEvents, const EventList& mcEvents,
                 const std::function<double( const Event& )>& fcn );
    Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
                   const std::function<double( const Event& )>& fcn, const unsigned int& minEvents = 10 );

    Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
                   const std::function<double( const Event& )>& fcn, const std::string& filename );
  };
} // namespace AmpGen

#endif /* end of include guard: AMPGEN_CHI2ESTIMATOR_H */
