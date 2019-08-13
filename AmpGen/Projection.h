#ifndef AMPGEN_PROJECTION_H
#define AMPGEN_PROJECTION_H

#include <stddef.h>
#include <functional>
#include <string>
#include <utility>

#include "TH1D.h"
#include "TH2D.h"
//#include "AmpGen/Event.h"

class TH1D;
class TH2D;

namespace AmpGen
{
  class Projection2D;
  class Event;
  class EventList; 

  class Projection
  {
    public:
      Projection();
      template <class FCN>
        Projection( const FCN& fcn, const std::string& name,
            const std::string& xAxisTitle, const size_t& nBins, const double& min, const double& max,
            const std::string& units = "" ) : Projection( std::function< double( const Event& )>( fcn ), name, xAxisTitle, nBins, min, max, units ) {}
      Projection( const std::function<double( const Event& )>& fcn, const std::string& name,
          const std::string& xAxisTitle, const size_t& nBins, const double& min, const double& max,
          const std::string& units = "" );
      const std::string name() const  ; 
      double operator()( const Event& evt ) const ;
      TH1D* operator()(const EventList& evt) const; 
      TH1D* plot(const std::string& prefix="") const;

      std::function<size_t( const Event& evt )> binFunctor() const;
      void setRange( const double& min, const double& max ){ m_min = (min); m_max = (max) ; }

      friend class Projection2D;
    private:
      std::function<double( const Event& )> m_func;
      std::string m_name;
      std::string m_xAxisTitle;
      std::string m_units;
      size_t m_nBins;
      double m_min;
      double m_max;
      double m_width;
  };

  class Projection2D
  {
    friend class Projection;
    Projection xAxis;
    Projection yAxis;

    public:
    Projection2D( const Projection& _xAxis, const Projection& _yAxis ) : xAxis( _xAxis ), yAxis( _yAxis ) {}

    TH2D* plot(const std::string& prefix="") const;

    std::pair<double, double> operator()( const Event& evt ) const;
  };

} // namespace AmpGen

#endif
