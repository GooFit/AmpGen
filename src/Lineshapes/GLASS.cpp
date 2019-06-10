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
  const auto props   = ParticlePropertiesList::get( particleName );
  const auto mass    = Parameter( particleName + "_mass", props->mass() );
  const auto width0  = Parameter( particleName + "_width", props->width() );
  const auto s0      = mass*mass;
  const auto a       = Parameter( particleName + "::GLASS::a", 2.07 );
  const auto r       = Parameter( particleName + "::GLASS::r", 3.32 );
  const auto phiR    = Parameter( particleName + "::GLASS::phiR", 0.00 );
  const auto phiF    = Parameter( particleName + "::GLASS::phiF", 0.00 );
  const auto R       = Parameter( particleName + "::GLASS::R", 1.00 );
  const auto F       = Parameter( particleName + "::GLASS::F", 1.00 );
  const auto q2      = Q2(s, s1, s2);
  const auto q20     = Q2(s0, s1, s2);
  const auto q       = safe_sqrt( q2 );
  const auto i       = Constant(0,1);
  const auto width   = width0 * mass * sqrt( q2 / ( q20 * s ) ); 
  const auto rbw     = s0 * width0 /(sqrt(q20) * (mass*mass - s - i*mass*width ));
  const auto den     = ( 1.0 + 0.5*r*a*q2 - i * a * q);
  const auto resTerm = R * exp( i*(phiR + 2*phiF) ) * rbw * ( 1.0 + 0.5*r*a*q2 + i * a * q );
  const auto nrTerm  = F * exp( i*phiF ) * ( a * q * cos(phiF) + ( 1.0 + 0.5 * r * a * q2 ) * sin(phiF)) * sqrt(s) / q;
  const auto total   = (nrTerm + resTerm)/den;

  ADD_DEBUG(width  , dbexpressions); 
  ADD_DEBUG(rbw    , dbexpressions);
  ADD_DEBUG(den    , dbexpressions); 
  ADD_DEBUG(resTerm/den, dbexpressions); 
  ADD_DEBUG(nrTerm/den , dbexpressions); 
  ADD_DEBUG(total  , dbexpressions); 
  return total;
}
