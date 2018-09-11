#ifndef AMPGEN_DALITZINTEGRATOR_H
#define AMPGEN_DALITZINTEGRATOR_H

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <string>
#include <utility>

#include "AmpGen/EventList.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Projection.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/WrappedMultiTF1.h"
#include "TF2.h"
#include "TRandom3.h"

namespace AmpGen
{
  /** @class DalitzIntegrator
   *  Class for doing 2D integrals using the Square Dalitz Plot (SQDP) method
   */

  class DalitzIntegrator
  {
  private:
    double m_min;
    double m_max;
    double m_mother;
    double m_massA;
    double m_massB;
    double m_massC;
    TRandom3* m_rand;

  public:
    typedef std::pair<double, double> sqCo;

    DalitzIntegrator( const double& mother, const double& m0, const double& m1, const double& m2 );
    void set( const double& mother, const double& m0, const double& m1, const double& m2 );

    sqCo getCoordinates( const Event& evt ) const;

    double sqDp1( const Event& evt ) const;
    double sqDp2( const Event& evt ) const;

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
             },
             0, 1, 0, 1, 0 );
      ROOT::Math::WrappedMultiTF1 wf1( f );
      ROOT::Math::AdaptiveIntegratorMultiDim ig;
      ig.SetFunction( wf1 );
      ig.SetRelTolerance( 0.000001 );
      double xmin[] = {0, 0};
      double xmax[] = {1, 1};
      double v = ig.Integral( xmin, xmax );
      return v / ( m_mother * m_mother );
    }

    template <class FCN>
    void debug( FCN fcn ) const
    {
      sqCo pos = {m_rand->Uniform(), m_rand->Uniform()};
      double event[12];
      for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
      setEvent( pos, event );
      for ( unsigned int i = 0; i < 12; ++i ) INFO( "Evt[" << i << "] = " << event[i] );
      fcn.debug( event );
      INFO( "Value = " << fcn( event ) );
    }
    void setRandom( TRandom3* rnd ) { m_rand = rnd; }

    double pCalc( const double& E, const double& m ) { return sqrt( E * E - m * m ); }

    void setEvent( const sqCo& x, double* event ) const;

    void setMin();

    void setMother( const double& m );

    double getMAB( sqCo coords ) const;

    double J( const sqCo& coords ) const;

    TH1D* makePlot( const std::function<double(const double*)>& fcn, const Projection& projection,
                    const std::string& name, unsigned int nSamples = 1000000 )
    {
      auto plot = projection.plot();
      double event[12];
      for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
      AmpGen::Event evtCache( 12 );
      for ( unsigned int i = 0; i < nSamples; ++i ) {
        sqCo pos = {m_rand->Uniform(), m_rand->Uniform()};
        setEvent( pos, event );
        evtCache.set( event );
        plot->Fill( projection( evtCache ), J( pos ) * fcn( event ) );
      }
      return plot;
    }

    TH2D* makePlot( const std::function<double(const double*)>& fcn, const Projection2D& projection,
                    const std::string& name, unsigned int nSamples = 1000000 )
    {
      auto plot = projection.plot();
      double event[12];
      for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
      AmpGen::Event evtCache( 12 );

      for ( unsigned int i = 0; i < nSamples; ++i ) {
        sqCo pos = {m_rand->Uniform(), m_rand->Uniform()};
        setEvent( pos, event );
        evtCache.set( event );
        auto obs_cos = projection( evtCache );
        plot->Fill( obs_cos.first, obs_cos.second, J( pos ) * fcn( event ) );
      }
      return plot;
    }
  };
} // namespace AmpGen

#endif
