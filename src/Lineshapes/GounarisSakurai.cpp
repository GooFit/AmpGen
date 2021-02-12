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
  Expression GS = kFactor( mass, width0, dbexpressions ) * FormFactor / (D-s);
  if(lineshapeModifier.find("Omega") != std::string::npos){
                // this implements rho/omega mixing from https://arxiv.org/pdf/hep-ex/0112031.pdf
                auto props_omega = ParticlePropertiesList::get("omega(782)0");
                Expression mass_omega = Parameter("omega(782)0_mass", props_omega->mass() );
                //Expression radius_omega = Parameter("omega(782)0_radius", props_omega->radius() );
                Expression width0_omega = Parameter("omega(782)0_width", props_omega->width() );
                
                const Expression kF_omega = 1.;//kFactor( mass_omega, width0_omega ) ;
                const Expression BW_omega = kF_omega / ( mass_omega * mass_omega - s - 1i*mass_omega * width0_omega );
                
                Expression delta_Amp = Parameter( particleName + "_deltaOmegaAmp", 0.00157 );
                Expression delta_Phase = Parameter( particleName + "_deltaOmegaPhase", 12.6 );
                Expression delta = delta_Amp * exp(1i* M_PI / 180 * delta_Phase);
                
                return GS * (1.+s/(mass_omega*mass_omega) * delta * BW_omega);
           }
return GS;
}
