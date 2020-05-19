#ifndef AMPGEN_PROJECTION_H
#define AMPGEN_PROJECTION_H

#include <stddef.h>
#include <functional>
#include <string>
#include <utility>

#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/Types.h"
#include "AmpGen/LiteSpan.h"

namespace AmpGen
{
  class Projection2D;
  class Event;
  class EventList; 

  class Projection
  {
    using keyedFunctors = KeyedFunctors<double, Event>; 
    public:
      Projection();
      template <class FCN>
      Projection( const FCN& fcn, const std::string& name,
          const std::string& xAxisTitle, const size_t& nBins, const double& min, const double& max,
          const std::string& units = "" ) : Projection( std::function< double( const Event& )>( fcn ), name, xAxisTitle, nBins, min, max, units ) {}
      Projection( const std::function<double( const Event& )>& fcn, const std::string& name,
          const std::string& xAxisTitle, const size_t& nBins, const double& min, const double& max,
          const std::string& units = "" );
      const std::string name() const; 
      template <class eventlist_type, class... ARGS> TH1D* operator()(const eventlist_type& evts, const ARGS... args) const 
      {
        return projInternal(evts, ArgumentPack(args...) ); 
      } 
      template <class eventlist_type, class... ARGS> std::tuple<std::vector<TH1D*>, THStack*> operator()(const eventlist_type& evts, 
          const keyedFunctors& weightFunction, const ARGS... args ) const 
      {
        return projInternal(evts, weightFunction, ArgumentPack(args...) );
      }

      double operator()( const Event& evt ) const;
      
      TH1D* plot(const std::string& prefix="") const;

      std::function<int( const Event& evt )> binFunctor() const;
      void setRange( const double& min, const double& max )
      { 
        m_min = min; 
        m_max = max; 
        m_width = (m_max-m_min)/double(m_nBins);
      }

      friend class Projection2D;
  ///  private:
      template <class eventlist_type> 
      TH1D* projInternal(const eventlist_type&, const ArgumentPack&) const; 
      template <class eventlist_type> 
      std::tuple<std::vector<TH1D*>, THStack*> projInternal(const eventlist_type&, const keyedFunctors&, const ArgumentPack&) const; 
      std::function<double( const Event& )> m_func;
      std::string m_name       = {""};
      std::string m_xAxisTitle = {""};
      std::string m_units      = {""};
      size_t m_nBins = {0};
      double m_min   = {0};
      double m_max   = {0};
      double m_width = {0};
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
  namespace PlotOptions {
    DECLARE_ARGUMENT(LineColor     , int);
    DECLARE_ARGUMENT(DrawStyle     , std::string);
    DECLARE_ARGUMENT(Selection     , std::function<bool( const Event& )>);
    DECLARE_ARGUMENT(Prefix        , std::string);
    DECLARE_ARGUMENT(Norm          , double);
    DECLARE_ARGUMENT(AddTo         , THStack*);
    DECLARE_ARGUMENT(AutoWrite     , bool);
  }
} // namespace AmpGen

#endif
