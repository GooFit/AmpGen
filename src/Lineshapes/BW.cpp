#include <string>
#include <complex>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Particle.h"

using namespace AmpGen;
using namespace AmpGen::fcn; 
using namespace std::complex_literals;

DEFINE_LINESHAPE( FormFactor )
{
  auto props                  = ParticlePropertiesList::get( particleName );
  Expression radius           = Parameter( particleName + "_radius", props->radius() );
  Expression mass             = Parameter( particleName + "_mass"  , props->mass() );
  const Expression q2         = make_cse( Q2( s, s1, s2 ) );
  const Expression q20        = make_cse( Q2( mass*mass, s1, s2 ) );
  int Lp = L;
  if( lineshapeModifier == "L0" ) Lp = 0;
  if( lineshapeModifier == "L1" ) Lp = 1;
  if( lineshapeModifier == "L2" ) Lp = 2;
  if( lineshapeModifier == "L3" ) Lp = 3;
  
  Expression                              FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, Lp ) );
  if ( lineshapeModifier == "BL" )        FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, Lp ) );
  if ( lineshapeModifier == "NFF" )       FormFactor = 1; 
  if ( lineshapeModifier == "BELLE2018" ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, Lp ) );

  if( L != 0 ){
    ADD_DEBUG( q2      , dbexpressions );
    ADD_DEBUG( radius  , dbexpressions );
  }
  ADD_DEBUG( FormFactor, dbexpressions );
  return FormFactor ; 
}

DEFINE_LINESHAPE( ExpFF )
{
  auto props                  = ParticlePropertiesList::get( particleName );
  Expression radius           = Parameter( particleName + "_radius", props->radius() );
  const Expression q2         = Q2( s, s1, s2 );
  const Expression FormFactor = Exp( -q2 * radius * radius / 2. );
  if( L != 0 ){
    ADD_DEBUG( radius, dbexpressions );
    ADD_DEBUG( s1, dbexpressions );
    ADD_DEBUG( s2, dbexpressions );
    ADD_DEBUG( s, dbexpressions );
    ADD_DEBUG( q2, dbexpressions );
  }
  ADD_DEBUG( FormFactor, dbexpressions ); 
  return FormFactor ;
}

DEFINE_LINESHAPE( None ) { return Constant( 1); }

DEFINE_LINESHAPE( BW )
{
  auto s_cse = make_cse(s);
  auto props = ParticlePropertiesList::get( particleName );
  const Expression& mass     = Parameter( particleName + "_mass", props->mass() );
  const Expression& width0   = Parameter( particleName + "_width", props->width() );
  const Expression& radius   = Parameter( particleName + "_radius", props->radius() );
  const Expression q2        = make_cse( Abs(Q2( s_cse, s1, s2 ) ) ) ;
  const Expression q20       = make_cse( Abs(Q2( mass * mass, s1, s2 )) );
  Expression                              FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( lineshapeModifier == "BL" )        FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  if ( lineshapeModifier == "BELLE2018" ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, L ) );
  Expression runningWidth                     = width( s_cse, s1, s2, mass, width0, radius, L, dbexpressions );
  const Expression BW = FormFactor / ( mass * mass - s_cse  -1i * mass * runningWidth );
  const Expression kf = kFactor( mass, width0, dbexpressions );
  ADD_DEBUG( FormFactor, dbexpressions );
  ADD_DEBUG( runningWidth, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  ADD_DEBUG( kf, dbexpressions );
  return lineshapeModifier == "BELLE2018" ? BW : kf*BW;
}

DEFINE_LINESHAPE( SBW )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression width0 = Parameter( particleName + "_width", props->width() );
  const Expression kF = kFactor( mass, width0 ) ;
  const Expression BW = 1 / ( mass * mass - s - 1i*mass * width0 );
  ADD_DEBUG( kF, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return kF * BW;
}
