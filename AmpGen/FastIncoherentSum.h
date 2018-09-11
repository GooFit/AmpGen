#ifndef AMPGEN_FASTINCOHERENTSUM_H
#define AMPGEN_FASTINCOHERENTSUM_H

#include <complex>
#include <iomanip>
#include <string>
#include <vector>

#include "AmpGen/EventList.h"
/// AmpGen
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"

namespace AmpGen
{
  class EventType;
  class FitFraction;
  class LinearErrorPropagator;
  class MinuitParameterSet;

  class FastIncoherentSum : public FastCoherentSum
  {

  public:
    FastIncoherentSum( const EventType& finalStates, AmpGen::MinuitParameterSet& mps,
                       const std::string& prefix = "Inco" );

    double getVal( const Event& evt ) const
    {
      double value( 0. );
      for ( auto& mE : m_matrixElements ) {
        value += std::norm( mE.coefficient * evt.getCache( mE.addressData ) );
      }
      return value;
    }
    double operator()( const Event& evt ) const { return prob( evt ); }
    double prob( const Event& evt ) const
    {
      DEBUG( "global weight = " << m_weight << ", pdf value = " << getVal( evt ) << ", norm = " << m_norm );
      return m_weight * getVal( evt ) / m_norm;
    }
    double prob_unnormalised( const Event& evt ) const { return getVal( evt ); }
    void prepare();
    double norm() const;
    double norm( const Bilinears& norms ) const;
    std::complex<double> norm( const unsigned int& i ) { return m_normalisations.get( i, 0 ); }
    //// warning - these functions are not properly tested ///
    std::vector<FitFraction> fitFractions( const LinearErrorPropagator& linProp );
  };
} // namespace AmpGen

#endif
