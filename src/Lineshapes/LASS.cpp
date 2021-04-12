#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Units.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

DEFINE_LINESHAPE( LASS )
{
  const auto props     = ParticlePropertiesList::get( particleName );
  Expression mass      = Parameter( particleName + "_mass", props->mass() );
  Expression width0    = Parameter( particleName + "_width", props->width() );
  Expression s0        = mass*mass;
  Expression a         = Parameter( "LASS::a", 2.07 );
  Expression r         = Parameter( "LASS::r", 3.32 );
  Expression q2        = Q2(s, s1, s2);
  Expression q20       = Q2(s0, s1, s2);
  Expression q         = safe_sqrt( q2 );
  Expression rho       = q/sqrt(s);
  Expression rho0      = sqrt(q20)/mass;
  Expression gamma     = width0 * rho / rho0;

  Expression i         = Constant(0,1);

  Expression tan_phi1   = 2*a*q/(2 + a*r*q2);
  Expression tan_phi2   = mass*gamma / (s0 - s);
  Expression phase      = atan(tan_phi1) + atan(tan_phi2);
  
  Expression bw           = ( mass * width0 ) / ( rho0 * ( s0 -s -i*mass*gamma ) );
  Expression nrPhaseShift = ( 2 + a * r * q2 + i * 2*a*q ) / ( 2+a*r*q2 - i*2*a*q );
  Expression NR           = 2*a*sqrt(s)/( 2 + a*r*q*q -2*i*a*q );
  ADD_DEBUG(s, dbexpressions );
  ADD_DEBUG(nrPhaseShift*bw, dbexpressions );
  ADD_DEBUG(NR, dbexpressions );
  ADD_DEBUG(NR + bw*nrPhaseShift, dbexpressions );
  ADD_DEBUG(sin(phase) * exp( i * phase ) , dbexpressions ); 
  ADD_DEBUG(sin(phase) * exp( i * phase ) / rho, dbexpressions );
  
  // Takes the splitting of the LASS into a "resonant" and "nonresonant"
  // component from Laura++, as a better implementation of the "gLASS" lineshape
  // the BW and NR amplitudes can be included separately with separate amplitudes and 

  if( lineshapeModifier == "BW" ) return bw * nrPhaseShift; 
  if( lineshapeModifier == "NR" ) return NR;
  return NR + bw*nrPhaseShift;

}
