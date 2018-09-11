#include <cmath>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"

//// implementation of Gounaris Sakurai Lineshape ///
using namespace AmpGen;

Expression safe_sqrt( const Expression& expression )
{
  return Ternary( expression > 0 , Sqrt(expression) , 0 );
}

Expression k( Expression mpipi )
{
  Constant mpi( 139.57018 );
  return safe_sqrt( 0.25 * mpipi - Pow( mpi, 2 ) );
}

Expression k2( Expression mpipi )
{
  Constant mpi( 139.57018 );
  return 0.25 * mpipi - Pow( mpi, 2 );
}

Expression kprime( Expression mpipi )
{
  Constant mpi( 139.57018 );
  return 1.0 / ( 8.0 * safe_sqrt( 0.25 * mpipi - Pow( mpi, 2 ) ) );
}

Expression h( Expression mpipi )
{
  Constant mpi( 139.57018 );
  Expression Logterm = Log( ( Sqrt( mpipi ) + 2.0 * k( mpipi ) ) / ( 2.0 * mpi ) );
  return ( 2.0 * k( mpipi ) * Logterm ) / ( M_PI * Sqrt( mpipi ) );
}
Expression hprime( Expression mpipi )
{
  return h( mpipi ) * ( 0.125 / ( k( mpipi ) * k( mpipi ) ) - 0.5 / mpipi ) + 0.5 / ( M_PI * mpipi );
}

Expression d( const Expression& mpipi, const Expression& m, DebugSymbols* dbexpressions = nullptr )
{
  Constant mpi( 139.57018 );
  Expression Logfactor = Log( ( m + 2 * k( m * m ) ) / ( 2.0 * mpi ) );
  Expression term1     = ( 3.0 / M_PI ) * Pow( mpi / k( m * m ), 2 ) * Logfactor;
  Expression term2     = m / ( 2.0 * M_PI * k( m * m ) );
  Expression term3     = Pow( mpi, 2 ) * m / ( M_PI * Pow( k( m * m ), 3 ) );

  return term1 + term2 - term3;
}
DEFINE_LINESHAPE( GounarisSakurai )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression radius = Parameter( particleName + "_radius", props->radius() );
  Expression width0 = Parameter( particleName + "_width", props->width() );

  Expression s0        = mass * mass;
  Expression GsTerm    = k2( s ) * ( h( s ) - h( s0 ) ) + k2( s0 ) * hprime( s0 ) * ( s0 - s );
  Expression BWReal    = s0 - s + width0 * s0 * GsTerm / ( k( s0 ) * k( s0 ) * k( s0 ) );
  Expression BWIm      = -mass * width( s, s1, s2, mass, width0, radius, L );
  const Expression q2  = Abs( Q2( s, s1, s2 ) );
  const Expression q20 = Abs( Q2( mass * mass, s1, s2 ) );
  const Expression J   = Constant(0,1);
  Expression BF = Sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  ADD_DEBUG( BF, dbexpressions );
  ADD_DEBUG( Sqrt( s ), dbexpressions );
  ADD_DEBUG( Sqrt( s0 ), dbexpressions );
  ADD_DEBUG( GsTerm, dbexpressions );
  ADD_DEBUG( BWReal, dbexpressions );
  ADD_DEBUG( BWIm, dbexpressions );
  ADD_DEBUG( h( s ), dbexpressions );
  ADD_DEBUG( h( s0 ), dbexpressions );
  ADD_DEBUG( k( s ), dbexpressions );
  ADD_DEBUG( k( s0 ), dbexpressions );
  ADD_DEBUG( d( s, mass ), dbexpressions );
  ADD_DEBUG( d( s0, mass ), dbexpressions );
  return  1000. * 1000. * kFactor( mass, width0, dbexpressions ) * BF / ( BWReal + J* BWIm );
}
