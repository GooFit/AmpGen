#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Spline.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

DEFINE_LINESHAPE( DecaySpline )
{
  const Expression q2 = Q2( s, s1, s2 ) / ( GeV * GeV );
  if ( dbexpressions != nullptr ) {
    ADD_DEBUG( q2, dbexpressions );
  }
  return getSpline( particleName, q2, "FF" );
}

DEFINE_LINESHAPE( GSpline )
{

  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() ) * GeV ;
  Expression radius = Parameter( particleName + "_radius", props->radius() ) * GeV ;
  Expression width0 = Parameter( particleName + "_width", props->width() ) * GeV;

  const Expression q2     = Abs( Q2( s, s1, s2 ) );
  const Expression sInGeV = s / ( GeV * GeV );

  bool useEFF            = lineshapeModifier.find( "EFF" ) != std::string::npos;
  bool useDispersiveTerm = lineshapeModifier.find( "Dis" ) != std::string::npos;
  bool useGFF            = lineshapeModifier.find( "GFF" ) != std::string::npos;

  Expression BF = fcn::sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) ) ;
  if( useEFF ) BF =  fcn::exp( -q2 * radius * radius / 2. );
  if ( useGFF ) {
    const Expression alpha = Parameter( particleName + "_alpha" );
    BF = fcn::exp( -( sInGeV - mass*mass ) * ( sInGeV - mass*mass )/( 2*alpha*alpha ) );
  }

  const Expression width_shape  = getSpline( particleName, sInGeV   , "Gamma", dbexpressions, true );
  const Expression width_norm   = getSpline( particleName, mass*mass, "Gamma", dbexpressions, true );
  const Expression norm         = kFactor( mass, width0 ) * BF;
  const Expression runningWidth = width0 * width_shape / width_norm;

  Expression real_part    = mass * mass;
  if( useDispersiveTerm ) real_part = real_part + getSpline( particleName, sInGeV, "dm2", dbexpressions, true ) - getSpline( particleName, mass*mass, "dm2", dbexpressions, true );
  const Expression self_energy = real_part - Constant(0,1) * mass * runningWidth; 

  const Expression BW        = norm / ( self_energy - s );

  ADD_DEBUG( mass, dbexpressions );
  ADD_DEBUG( q2, dbexpressions );
  ADD_DEBUG( BF, dbexpressions );
  ADD_DEBUG( s, dbexpressions );
  ADD_DEBUG( width_shape, dbexpressions );
  ADD_DEBUG( width_norm, dbexpressions );
  ADD_DEBUG( runningWidth, dbexpressions );
  ADD_DEBUG( kFactor( mass, width0 ), dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return BW;
}

DEFINE_LINESHAPE( FormFactorSpline )
{
  auto props                    = ParticlePropertiesList::get( particleName );
  Expression mass               = Parameter( particleName + "_mass", props->mass() );
  Expression radius             = Parameter( particleName + "_radius", props->radius() );
  Expression width0             = Parameter( particleName + "_width", props->width() );
  const Expression q2           = Abs( Q2( s, s1, s2 ) );
  const Expression q20          = Abs( Q2( mass * mass, s1, s2 ) );
  const Expression runningWidth = width( s, s1, s2, mass, width0, radius, L, dbexpressions );
  const Expression BW =
    GeV * getSpline( particleName, s, "FF" ) / ( -s + mass * mass  - Constant(0,1) * mass * runningWidth );
  return BW;
}

DEFINE_LINESHAPE( MIPWA )
{
  bool cont       = lineshapeModifier.find( "continue" ) != std::string::npos;
  Expression real = getSpline( particleName, s / ( GeV * GeV ), "Re", dbexpressions, cont );
  Expression imag = getSpline( particleName, s / ( GeV * GeV ), "Im", dbexpressions, cont );
  ADD_DEBUG( real, dbexpressions );
  ADD_DEBUG( imag, dbexpressions );
  ADD_DEBUG( s, dbexpressions );
  Expression J = Constant(0,1);
  return lineshapeModifier != "Polar" ? ( real + J*imag ) : real * ( Cos( imag ) + J*Sin( imag ) );
}

DEFINE_LINESHAPE( InelasticSpline )
{
  Expression inelasticity = getSpline( particleName, s / ( GeV* GeV), "Eta", dbexpressions, true );
  Expression phaseShift   = getSpline( particleName, s / ( GeV* GeV), "Delta", dbexpressions, true );
  Expression J = Constant(0,1);
  Expression phi             = Cos( 2 * phaseShift ) + J * Sin( 2 * phaseShift ) ;
  return ( inelasticity * phi - 1 ) * ( -0.5*J ) ;
}
