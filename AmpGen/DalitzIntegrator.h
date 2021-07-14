#ifndef AMPGEN_DALITZINTEGRATOR_H
#define AMPGEN_DALITZINTEGRATOR_H

#include <stddef.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <string>
#include <utility>

#include "AmpGen/NumericalIntegration.h"

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
      typedef std::pair<double, double> sqCo;

      DalitzIntegrator( const double& s0, const double& s1, const double& s2, const double& s3);

      template <typename FCN> double integrateDP( FCN fcn, const double& s) const;

      template <typename FCN> double integrateDP( FCN fcn ) const
      {
        return integrateDP( fcn, m_s0 );
      }
      double getMAB(sqCo coords)   const;
      double J(const sqCo& coords) const;
      double getMAB(sqCo coords  , const double& s) const;
      double J(const sqCo& coords, const double& s) const;
      double sqDp1(const Event& evt) const;
      double sqDp2(const Event& evt) const;
      double safe_sqrt(const double& x ) const { return x > 0 ? sqrt(x) : 0; }
      void setEvent(const sqCo& x, double* event, const double& s) const;
      void debug() const; 
      void setEvent(const sqCo& x, double* event) const;
      void set(const double& s0, const double& s1, const double& s2, const double& s3);
      void setMin();
      void setMother(const double& s);

      TH1D* makePlot( const std::function<double(const double*)>& fcn, const Projection& projection,
          const std::string& name, const size_t& nSamples = 1000000 );

      TH2D* makePlot( const std::function<double(const double*)>& fcn, const Projection2D& projection,
          const std::string& name, const size_t& nSamples = 1000000 );
      
      sqCo getCoordinates( const Event& evt ) const;

//    private:
      double    m_min;
      double    m_max;
      double    m_s0;
      double    m_s1;
      double    m_s2;
      double    m_s3;
   //   double integrate_internal( TF2& fcn ) const ;
   
  };
  template <typename FCN>
  double DalitzIntegrator::integrateDP( FCN fcn, const double& s) const
  {
    double event[12];
    for(size_t i = 0; i < 12; ++i) event[i] = 0; 
    return integrate<2>(  [&]( const std::array<double,2>& x ) {
      setEvent( *reinterpret_cast<const sqCo*>(&x) , event, s);
      return J( *reinterpret_cast<const sqCo*>(&x) , s) * std::real( fcn(event) ); }, std::array<double,2>{0., 0.}, std::array<double, 2>{1.,1.} ) /s;
  }

} // namespace AmpGen

#endif
