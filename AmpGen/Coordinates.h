#ifndef AMPGEN_COORDINATES_H
#define AMPGEN_COORDINATES_H 1

#include "AmpGen/EventList.h"
#include "AmpGen/Kinematics.h"

#include <functional>
#include <vector>

namespace AmpGen
{
  namespace functors
  {
    std::function<std::vector<double>( const AmpGen::Event& )> makeOptimisedFunctors( const double& units = 1000 )
    {
      auto kstar_hcos = HelicityCosine( {0}, {0, 1, 2, 3}, {0, 1} ); /// kstarHelicityCosine  ////
      auto rho_hcos   = HelicityCosine( {2}, {0, 1, 2, 3}, {2, 3} ); /// rhoHelicityCosine ////
      return [units, kstar_hcos, rho_hcos]( const Event& event ) -> std::vector<double> {
        return {PHI( event ), sqrt( event.s( 0, 1 ) ) / units, sqrt( event.s( 2, 3 ) ) / units, kstar_hcos( event ),
                rho_hcos( event )};
      };
    }

    std::function<std::vector<double>( const Event& )> makeOptimisedFunctorsSym( const double& units = 1000 )
    {
      auto kstar_hcos = HelicityCosine( {0}, {0, 1, 2, 3}, {0, 1} ); /// kstarHelicityCosine  ////
      auto rho_hcos   = HelicityCosine( {2}, {0, 1, 2, 3}, {2, 3} ); /// rhoHelicityCosine ////

      auto kstar_hcos_prime = HelicityCosine( {0}, {0, 2, 1, 3}, {0, 2} ); /// kstarHelicityCosine  ////
      auto rho_hcos_prime   = HelicityCosine( {1}, {0, 2, 1, 3}, {1, 3} ); /// rhoHelicityCosine ////

      return
          [units, kstar_hcos, rho_hcos, kstar_hcos_prime, rho_hcos_prime]( const Event& event ) -> std::vector<double> {
            double s01 = sqrt( event.s( 0, 1 ) );
            double s02 = sqrt( event.s( 0, 2 ) );
            if ( s01 > s02 )
              return {phi( event, 0, 1, 3, 2 ), sqrt( event.s( 0, 2 ) ) / units, sqrt( event.s( 1, 3 ) ) / units,
                      kstar_hcos_prime( event ), rho_hcos_prime( event )};
            else
              return {phi( event, 0, 2, 3, 1 ), sqrt( event.s( 0, 1 ) ) / units, sqrt( event.s( 2, 3 ) ) / units,
                      kstar_hcos( event ), rho_hcos( event )};
          };
    }

    std::function<bool( const Event& )> makeKsVeto( const double& width = 0.015, const double& units = 1000 )
    {
      return [width, units]( const Event& evt ) -> bool {
        double m1 = sqrt( evt.s( 2, 3 ) );
        double m2 = sqrt( evt.s( 1, 3 ) );
        return ( abs( m1 - 0.4976 * units ) < width || abs( m2 - 0.4976 * units ) < width );
      };
    }
  } // namespace functors
} // namespace AmpGen

#endif
