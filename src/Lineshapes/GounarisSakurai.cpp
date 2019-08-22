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
using namespace std::complex_literals;

Expression q2( const Expression& s )
{
  double mpi = ParticlePropertiesList::get("pi+")->mass();
  return s - 4*mpi*mpi;
}

Expression q( const Expression& s )
{
  return safe_sqrt(q2(s));
}

Expression logTerm( const Expression& s )
{
  double mpi = ParticlePropertiesList::get("pi+")->mass();
  return log( ( sqrt(s) + safe_sqrt(q2(s)) ) / ( 2*mpi ) );
}

DEFINE_LINESHAPE( GounarisSakurai )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression radius = Parameter( particleName + "_radius", props->radius() );
  Expression width0 = Parameter( particleName + "_width", props->width() );

  Expression s0        = mass * mass;
   
  Expression t1 = 2 * s0 * q2(s) * q(s) *logTerm(s)/(sqrt(s)*fpow(q(s0),3));
  Expression t2 = logTerm(s0)*( q2(s)*(s - 3*s0) + s*(s0-s) )/(mass*q2(s0)) ;
  Expression t3 = (s0-s)/q(s0);
  
  Expression M2 = s0 + width0 * (t1 + t2 + t3 )/M_PI; 
  bool useRadius    = lineshapeModifier.find("IncludeRadiusInWidth") != std::string::npos;
  bool normalisedFF = lineshapeModifier.find("BL")                   != std::string::npos;
  Expression  D = M2 - 1i*mass*width(s, s1, s2, mass, width0, useRadius ? radius : 0, L);
  const Expression q2  = abs( Q2( s, s1, s2 ) );
  Expression FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( normalisedFF ) FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  ADD_DEBUG( FormFactor   , dbexpressions);
  ADD_DEBUG( M2 , dbexpressions);
  ADD_DEBUG( D, dbexpressions);
  return kFactor( mass, width0, dbexpressions ) * FormFactor / (D-s);
}
