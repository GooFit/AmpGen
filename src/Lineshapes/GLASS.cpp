#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Units.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

DEFINE_LINESHAPE( GLASS )
{
  WARNING("GLASS lineshape does not work as expected");
  const auto props     = ParticlePropertiesList::get( particleName );
  Expression mass      = Parameter( particleName + "_mass", props->mass() );
  Expression width0    = Parameter( particleName + "_width", props->width() );
  Expression s0        = mass*mass;
  Expression a         = Parameter( particleName + "::GLASS::a", 2.07 );
  Expression r         = Parameter( particleName + "::GLASS::r", 3.32 );
  Expression phi_R     = Parameter( particleName + "::GLASS::phiR", 0.00 );
  Expression phi_F     = Parameter( particleName + "::GLASS::phiF", 0.00 );
  Expression R         = Parameter( particleName + "::GLASS::R", 1.00 );
  Expression F         = Parameter( particleName + "::GLASS::F", 1.00 );
  Expression q2        = Q2(s, s1, s2);
  Expression q20       = Q2(s0, s1, s2);
  Expression q         = safe_sqrt( q2 );
  Expression rho       = q/sqrt(s);
  Expression rho0      = sqrt(q20)/mass;
  Expression gamma     = width0 * rho / rho0;
  Expression dF        = make_cse(M_PI*phi_F/180. + atan(2*a*q/(2 + a*r*q2)));
  Expression dR        = make_cse(M_PI*phi_R/180. + atan(mass*gamma / (s0 - s)));
  Expression i         = Constant(0,1);
  ADD_DEBUG(s, dbexpressions );
  ADD_DEBUG(s0, dbexpressions );
  ADD_DEBUG(dF, dbexpressions ); 
  ADD_DEBUG(dR, dbexpressions );
  ADD_DEBUG(phi_F, dbexpressions ); 
  ADD_DEBUG(phi_R, dbexpressions ); 
  ADD_DEBUG(gamma, dbexpressions ); 
  
  return R*sin(dR)*exp(i*dR)*exp(2*i*dF) + F*sin(dF)*exp(i*dF);
}
