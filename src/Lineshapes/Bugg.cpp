#include <cmath>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"

//// implementation of scalar S-wave Lineshape suggested in D. V. Bugg, J. Phys. G34:151, 2007, arXiv:hep-ph/0608081 ///

using namespace AmpGen;

Expression rho( const Expression& s, const Expression& s0 )
{
  Expression a2 = 1.0 - ( 4 * s0 ) / s;
  return Ternary( a2 > 0, a2, Constant( 0 ) );
}

Expression rho_2( const Expression& s, const Expression& s0 )
{
  return fcn::complex_sqrt( 1.0 - ( 4 * s0 ) / s );
}

Expression q( const Expression& s, const Expression& s0 ) { return Ternary( s > 4 * s0, s - 4 * s0, 4 * s0 - s ); }

Expression rho_4pi( const Expression& s, const Expression& lambda, const Expression& s0 )
{
  return 1. / ( 1 + Exp( lambda * ( s0 - s ) ) );
}

Expression Buggj1( const Expression& s, const Expression& m0 )
{

  Expression rho_pipi = Sqrt( rho( s, m0 * m0 ) );
  Expression addTerm  = Ternary( rho_pipi > 0., rho_pipi * Log( ( 1. - rho_pipi ) / ( 1. + rho_pipi ) ), 0. );
  return ( 2. + addTerm ) / M_PI;
}

Expression Gamma_4pi( const Expression& s, const Expression& m0, const Expression& M, const Expression& g_4pi,
                      const Expression& lambda, const Expression& s0 )
{
  return Ternary( s > 16. * m0 * m0, g_4pi * rho_4pi( s, lambda, s0 ) / rho_4pi( M * M, lambda, s0 ), 0. );
}

DEFINE_LINESHAPE( Bugg )
{
  Expression M      = Parameter( "Bugg::M", 0.935 );
  Expression b1     = Parameter( "Bugg::b1", 1.302 );
  Expression b2     = Parameter( "Bugg::b2", 0.340 );
  Expression A      = Parameter( "Bugg::A", 2.426 );
  Expression g_4pi  = Parameter( "Bugg::g_4pi", 0.011 );
  Expression g_2K   = Parameter( "Bugg::g_2K", 0.6 );
  Expression g_2eta = Parameter( "Bugg::g_2eta", 0.2 );
  Expression alpha  = Parameter( "Bugg::alpha", 1.3 );
  Constant mPiPlus( 0.139570 );
  Constant mKPlus( 0.493677  );
  Constant mEta( 0.547863  );
  Expression J = Constant(0,1);

  Expression sA         = Parameter( "Bugg::sA", 0.41 ) * mPiPlus * mPiPlus;
  Expression s0_4pi     = Parameter( "Bugg::s0_4pi", 7.082 / 2.845 );
  Expression lambda_4pi = Parameter( "Bugg::lambda_4pi", 2.845 );

  Expression z = Buggj1(s, mPiPlus) - Buggj1(M * M, mPiPlus);

  Expression g1sg      = M * ( b1 + b2 * s ) * Exp( -( s - M * M ) / A );
  Expression adlerZero = ( s - sA ) / ( M * M - sA );

  Expression gamma_2pi = g1sg * adlerZero * rho_2( s, mPiPlus * mPiPlus );

  Expression gamma_2K = g_2K * g1sg * s / ( M * M ) * Exp( -alpha * q( s, mKPlus * mKPlus ) ) *
                        rho_2( s, mKPlus * mKPlus );

  Expression gamma_2eta =
      g_2eta * g1sg * s / ( M * M ) * Exp( -alpha * q( s, mEta * mEta ) ) * rho_2( s, mEta * mEta );
  Expression gamma_4pi = M * Gamma_4pi( s, mPiPlus, M, g_4pi, lambda_4pi, s0_4pi );

  Expression Gamma_tot = gamma_2pi + gamma_2K + gamma_2eta + gamma_4pi;
  Expression az        = g1sg * adlerZero;

  Expression iBW = M * M - s - az * z - J * Gamma_tot;

  Expression BW = 1. / iBW;
  ADD_DEBUG( M, dbexpressions );
  ADD_DEBUG( g1sg, dbexpressions );
  ADD_DEBUG( gamma_2pi, dbexpressions );
  ADD_DEBUG( gamma_2K, dbexpressions );
  ADD_DEBUG( gamma_2eta, dbexpressions );
  ADD_DEBUG( gamma_4pi, dbexpressions );
  ADD_DEBUG( sA, dbexpressions );
  ADD_DEBUG( adlerZero, dbexpressions );
  ADD_DEBUG( Gamma_tot, dbexpressions );
  ADD_DEBUG( s, dbexpressions );
  ADD_DEBUG( z, dbexpressions );
  ADD_DEBUG( Buggj1( s, mPiPlus ), dbexpressions );
  ADD_DEBUG( iBW, dbexpressions );
  ADD_DEBUG( az * z, dbexpressions );
  ADD_DEBUG( rho_4pi( s, lambda_4pi, s0_4pi ), dbexpressions ); 
  return lineshapeModifier == "Az" ? BW * ( s - sA ) * g1sg : BW;
}
