#ifndef AMPGEN_DALITZINTEGRATOR_H
#define AMPGEN_DALITZINTEGRATOR_H

#include <stddef.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <string>
#include <utility>

#include "AmpGen/simd/utils.h"
#include "AmpGen/NumericalIntegration.h"
#include "AmpGen/ProfileClock.h"
class TH1D;
class TH2D;

namespace AmpGen
{
  class Event;
  class Projection2D;
  class Projection;
  /** @class DalitzIntegrator
   *  Class for doing 2D integrals using the Square Dalitz Plot (SQDP) method
   */

  class DalitzIntegrator
  {
    public:
      typedef std::pair<float_v, float_v> sqCo;

      DalitzIntegrator( const double& s0, const double& s1, const double& s2, const double& s3);

      template <typename FCN> double integrateDP( FCN&& fcn, const double& s) const;
      template <typename FCN> double integrateDP( FCN&& fcn ) const
      {
        return integrateDP( fcn, m_s0 );
      }
      float_v getMAB(const sqCo& coords)   const;
      float_v J(const sqCo& coords) const;
      float_v getMAB(const sqCo& coords  , const double& s) const;
      float_v J(const sqCo& coords, const double& s) const;
      double sqDp1(const Event& evt) const;
      double sqDp2(const Event& evt) const;
      float_v safe_sqrt(const float_v& x ) const { 
        #if ENABLE_AVX 
        return select( x > 0. , sqrt(x) , 0. ); 
        #else 
          return x > 0 ? sqrt(x) : 0;
        #endif
      }
      void setEvent(const sqCo& x, float_v* event, const double& s) const;
      void debug() const; 
      void setEvent(const sqCo& x, float_v* event) const;
      void set(const double& s0, const double& s1, const double& s2, const double& s3);
      void setMin();
      void setMother(const double& s);

      TH1D* makePlot( const std::function<double(const double*)>& fcn, const Projection& projection,
          const std::string& name, const size_t& nSamples = 1000000 );

      TH2D* makePlot( const std::function<double(const double*)>& fcn, const Projection2D& projection,
          const std::string& name, const size_t& nSamples = 1000000 );
      
      sqCo getCoordinates( const Event& evt ) const;

    private:
      double    m_min;
      double    m_max;
      double    m_s0;
      double    m_s1;
      double    m_s2;
      double    m_s3;
   
  };
  template <typename FCN>
  double DalitzIntegrator::integrateDP( FCN&& fcn, const double& s) const
  {
    float_v event[12] = {0.};
    ProfileClock pc1; 
    auto i1 = integrate<2>(  [&]( const std::array<float_v,2>& x ) {
      setEvent( *reinterpret_cast<const sqCo*>(&x) , event, s);
      return J( *reinterpret_cast<const sqCo*>(&x) , s) * real( fcn(event) ); }, std::array<double,2>{0., 0.}, std::array<double, 2>{1.,1.} ) /s;
    return i1; 
  }


} // namespace AmpGen

#endif
