#include <math.h>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"

//// implementation of Gounaris Sakurai Lineshape ///
using namespace AmpGen;
using namespace AmpGen::fcn;

Expression q2( const Expression& s )
{
  double mpi = ParticlePropertiesList::get("pi+")->mass();
  return s - 4*mpi*mpi;
}

Expression q( const Expression& s )
{
  return safe_sqrt(q2(s));
}

Expression L( const Expression& s )
{
  double mpi = ParticlePropertiesList::get("pi+")->mass();
  return log( ( sqrt(s) + safe_sqrt(q2(s)) ) / ( 2*mpi ) );
}

Expression h( const Expression& s )
{
  return  q(s) * L(s) / (M_PI * sqrt(s) );
}
Expression hprime( const Expression& s )
{
  double mpi = ParticlePropertiesList::get("pi+")->mass();
  return ( L(s) * 4 * mpi * mpi /(q(s) * sqrt(s)) + 1) / (2*M_PI*s);
}
Expression M( const Expression& s, const Expression& width, const Expression& s0, const Expression& pion_mass ) {
  Expression pf = 2.*width*s0/pow(q(s0),3);
  Expression t1 = q2(s)  * ( h(s) - h(s0));
  Expression t2 = q2(s0) * ( L(s0) * 4 * pion_mass * pion_mass/(q(s0) * sqrt(s0)) + 1) ;
  return s0 + pf*(t1+t2*(s0-s)/(2.*M_PI*s0)); 
}

DEFINE_LINESHAPE( GounarisSakurai )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression radius = Parameter( particleName + "_radius", props->radius() );
  Expression width0 = Parameter( particleName + "_width", props->width() );

  Expression s0        = mass * mass;
    
  Expression BWIm      = -mass * width( s, s1, s2, mass, width0, 0, L );
  Expression BWRe      = M(s, width0, s0, ParticlePropertiesList::get("pi+")->mass() ) - s;
  const Expression q2  = abs( Q2( s, s1, s2 ) );
  const Expression J   = Constant(0,1);
  Expression FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( lineshapeModifier == "BL" ) FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  ADD_DEBUG( FormFactor   , dbexpressions);
  ADD_DEBUG( BWRe , dbexpressions);
  ADD_DEBUG( BWIm , dbexpressions);
  ADD_DEBUG( h(s) , dbexpressions);
  ADD_DEBUG( h(s0), dbexpressions);
  return kFactor( mass, width0, dbexpressions ) * FormFactor / ( BWRe + J* BWIm );
}
