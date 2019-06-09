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
  const auto props  = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression width0 = Parameter( particleName + "_width", props->width() );
  Expression s0     = mass*mass;
  Expression a      = Parameter( particleName + "::GLASS::a", 2.07 );
  Expression r      = Parameter( particleName + "::GLASS::r", 3.32 );
  Expression phiR   = Parameter( particleName + "::GLASS::phiR", 0.00 );
  Expression phiF   = Parameter( particleName + "::GLASS::phiF", 0.00 );
  Expression R      = Parameter( particleName + "::GLASS::R", 1.00 );
  Expression F      = Parameter( particleName + "::GLASS::F", 1.00 );
  Expression q2     = Q2(s, s1, s2);
  Expression q20    = Q2(s0, s1, s2);
  Expression q      = safe_sqrt( q2 );
  Expression i      = Constant(0,1);
  Expression mR     = mass; 
  Expression pR     = sqrt( q20 );
  Expression gammaR = width0;

  auto g = width0 * mass * sqrt( q2 / ( q20 * s ) ); 
  auto propagator_relativistic_BreitWigner = 1./(mass*mass - s - i*mass*g);
  auto cot_deltaF  = 1.0/(a*q) + 0.5*r*q;
  auto qcot_deltaF = 1.0/a + 0.5*r*q2;
  auto expi2deltaF = (qcot_deltaF + i * q) / (qcot_deltaF - i * q);
  auto resonant_term_T     = R * exp( i*(phiR + 2*phiF) ) * propagator_relativistic_BreitWigner * mR * gammaR * mR / pR * expi2deltaF;
  auto non_resonant_term_F = F * exp( i*phiF ) * (cos(phiF) + cot_deltaF * sin(phiF)) * sqrt(s) / (qcot_deltaF -i * q);
  auto LASS_contribution   = non_resonant_term_F + resonant_term_T;

  ADD_DEBUG( g                                  , dbexpressions); 
  ADD_DEBUG( propagator_relativistic_BreitWigner, dbexpressions); 
  ADD_DEBUG( cot_deltaF                         , dbexpressions); 
  ADD_DEBUG( qcot_deltaF                        , dbexpressions); 
  ADD_DEBUG( expi2deltaF                        , dbexpressions); 
  ADD_DEBUG( resonant_term_T                    , dbexpressions); 
  ADD_DEBUG( non_resonant_term_F                , dbexpressions); 
  ADD_DEBUG( LASS_contribution                  , dbexpressions); 
  return LASS_contribution;
}
