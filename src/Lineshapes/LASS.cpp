#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

Expression phase( const Expression& phi )
{ 
  return Cos(phi) + Constant(0,1)*Sin(phi);
}


DEFINE_LINESHAPE( LASS )
{
  const auto props                 = ParticlePropertiesList::get( particleName );
  const Expression mass            = Parameter( particleName + "_mass", props->mass() );
  const Expression radius          = Parameter( particleName + "_radius", props->radius() );
  const Expression width0          = Parameter( particleName + "_width", props->width() );
  const Expression a               = Parameter( "LASS::a", 2.07 );
  const Expression r               = Parameter( "LASS::r", 3.32 );
  const Expression g0              = Parameter( "LASS::g" + lineshapeModifier, 0.0 );
  const Expression q2              = Q2( s, s1, s2 ) / ( GeV*GeV );
  const Expression q20             = Q2( mass*mass, s1, s2 ) / (GeV*GeV) ;
  const Expression q               = fcn::safe_sqrt( q2 );
  const Expression gamma           = width( s, s1, s2, mass, width0, radius, 0, dbexpressions );
  const Expression scatteringPhase = ATan( 2 * a * q / ( 2 + a * r * q2 ) );
  const Expression resonancePhase  = ATan( mass * gamma / ( mass * mass - s ) );
  const Expression tphase          = scatteringPhase + resonancePhase;
  const Expression rho             = GeV * fcn::safe_sqrt( q2 / s );
  const Expression I               = Constant(0,1);
  const Expression returnValue     = Sin( tphase ) * phase(tphase) / rho;

  const Expression phi1            = ( 2 * a * q / ( 2 + a * r * q2 ) );
  const Expression phi2            = ( mass * gamma / ( mass * mass - s ) );

  const Expression rho0            = GeV * fcn::safe_sqrt( q20 / (mass*mass) );
  const Expression ow              = width0 * rho / rho0; 
  const Expression mG              = mass*gamma;
  const Expression mG0             = mass*width0/rho0;
  const Expression del             = mass*mass - s;

  const Expression phi1ByRho           = 2 * a * fcn::sqrt(s) / ( GeV* ( 2 + a * r *q2 ) );
  const Expression numerator       = mG0 + phi1ByRho*del;

  const Expression denom           = del - mG*phi1 - I*( mG + phi1*del );
  ADD_DEBUG( s, dbexpressions );
  ADD_DEBUG( q2, dbexpressions );
  ADD_DEBUG( q, dbexpressions );
  ADD_DEBUG( scatteringPhase, dbexpressions );
  ADD_DEBUG( resonancePhase, dbexpressions );
  ADD_DEBUG( rho, dbexpressions );
  ADD_DEBUG( returnValue, dbexpressions );
  ADD_DEBUG( gamma, dbexpressions );
  ADD_DEBUG( numerator / (denom), dbexpressions );


  Expression production_form_factor = 1;

  if ( lineshapeModifier == "pFF" ) {
    production_form_factor = 1 - Parameter( "LASS::pFF::f" ) * Exp( -q2 * Parameter( "LASS::pFF::alpha" ) );
  }
  ADD_DEBUG( production_form_factor, dbexpressions );
  return production_form_factor * numerator / denom;
}

DEFINE_LINESHAPE( gLASS )
{
  auto props           = ParticlePropertiesList::get( particleName );
  Expression mass      = Parameter( particleName + "_mass", props->mass() );
  Expression radius    = Parameter( particleName + "_radius", props->radius() );
  Expression width0    = Parameter( particleName + "_width", props->width() );
  const Expression sK  = 493.677 * 493.677;
  const Expression sPi = 139.57 * 139.57;
  Expression a         = Parameter( "gLASS::a" + lineshapeModifier, 2.07 );   /// short distance coupling
  Expression r         = Parameter( "gLASS::r" + lineshapeModifier, 3.32 );   /// range parameter
  Expression F         = Parameter( "gLASS::F" + lineshapeModifier, 1.0 );    /// vec[0];
  Expression R         = Parameter( "gLASS::R" + lineshapeModifier, 1.0 );    /// resonant coupling
  Expression phiF      = Parameter( "gLASS::phiF" + lineshapeModifier, 0.0 ); /// F-phase
  Expression phiS      = Parameter( "gLASS::phiS" + lineshapeModifier, 0.0 ); /// R-phase

  const Expression q2              = Q2( s, s1, s2 ) / ( GeV * GeV );
  const Expression q               = Sqrt( q2 );
  const Expression gamma           = width( s, s1, s2, mass, width0, radius, 0 );
  const Expression scatteringPhase = phiF + ATan( 2 * a * q / ( 2 + a * r * q2 ) );
  const Expression resonancePhase  = phiS + ATan( mass * gamma / ( mass * mass - s ) );
  const Expression rho             = GeV * Sqrt( q2 / s );
  const Expression returnValue     = ( F * Sin( scatteringPhase ) * phase(scatteringPhase)  +
      R * Sin( resonancePhase ) * phase( resonancePhase + 2 * scatteringPhase ) ) / rho;
  ADD_DEBUG( q, dbexpressions );
  ADD_DEBUG( rho, dbexpressions );
  ADD_DEBUG( gamma, dbexpressions );
  ADD_DEBUG( rho, dbexpressions );
  ADD_DEBUG( scatteringPhase, dbexpressions );
  ADD_DEBUG( resonancePhase, dbexpressions );

  return returnValue;
}
