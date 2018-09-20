#ifndef AMPGEN_DALITZINTEGRATOR_H
#define AMPGEN_DALITZINTEGRATOR_H
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <string>
#include <utility>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TRandom3.h>

#include "AmpGen/EventList.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Projection.h"

namespace AmpGen
{
  /** @class DalitzIntegrator
   *  Class for doing 2D integrals using the Square Dalitz Plot (SQDP) method
   */

  class DalitzIntegrator
  {
    private:
      double    m_min;
      double    m_max;
      double    m_s0;
      double    m_s1;
      double    m_s2;
      double    m_s3;
      double integrate_internal( TF2& fcn ) const ;
   
    public:
      typedef std::pair<double, double> sqCo;

      DalitzIntegrator( const double& s0, const double& s1, const double& s2, const double& s3);

      sqCo getCoordinates( const Event& evt ) const;

      template <class FCN>
      double integrate( FCN fcn ) const
      {
        double event[12];
        for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
        TF2 f( "fcn",
            [&]( double* x, double* p ) {
            sqCo pos = {x[0], x[1]};
            setEvent( pos, event );
            return J( pos ) * std::real( fcn( event ) );
            }, 0, 1, 0, 1, 0 );
        return integrate_internal(f) / m_s0;
      }
      
      template <class FCN>
      double integrate( FCN fcn, const double& s) const
      {
        double event[12];
        for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
        TF2 f( "fcn",
            [&]( double* x, double* p ) {
            sqCo pos = {x[0], x[1]};
            setEvent( pos, event, s);
            return J( pos, s ) * std::real( fcn( event ) );
            }, 0, 1, 0, 1, 0 );
        return integrate_internal(f) / s;
      }

      double getMAB( sqCo coords )   const;
      double J( const sqCo& coords ) const;
     
      double getMAB( sqCo coords  , const double& s )   const;
      double J( const sqCo& coords, const double& s ) const;
      void setEvent( const sqCo& x, double* event, const double& s ) const;

      double sqDp1( const Event& evt ) const;
      double sqDp2( const Event& evt ) const;
      double safe_sqrt( const double& x ) const { return x > 0 ? sqrt(x) : 0; }

      void debug() const; 
      void setEvent( const sqCo& x, double* event ) const;
      void set(         const double& s0, const double& s1, const double& s2, const double& s3);
      void setMin();
      void setMother( const double& s );

      TH1D* makePlot( const std::function<double(const double*)>& fcn, const Projection& projection,
          const std::string& name, const size_t& nSamples = 1000000 );

      TH2D* makePlot( const std::function<double(const double*)>& fcn, const Projection2D& projection,
          const std::string& name, const size_t& nSamples = 1000000 );
  };
} // namespace AmpGen

#endif
