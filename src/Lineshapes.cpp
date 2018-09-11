#include "AmpGen/Lineshapes.h"

#include <cmath>
#include <ostream>

#include "AmpGen/MsgService.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

Expression AmpGen::Q2( const Expression& Msq, const Expression& M1sq, const Expression& M2sq )
{
  const Expression num = Msq / 4. - ( M1sq + M2sq ) / 2. + ( M1sq - M2sq ) * ( M1sq - M2sq ) / ( 4 * Msq );
  return num;
}

Expression AmpGen::kFactor( const Expression& mass, const Expression& width, DebugSymbols* dbexpressions )
{
  const Expression massInGeV  = mass / GeV;
  const Expression widthInGeV = width / GeV;
  const Expression gammaInGeV = massInGeV * Sqrt( massInGeV * massInGeV + widthInGeV * widthInGeV );
  const Expression kInGeV =
      2 * sqrt( 2 ) * massInGeV * widthInGeV * gammaInGeV / ( M_PI * Sqrt( massInGeV * massInGeV + gammaInGeV ) );
  ADD_DEBUG( gammaInGeV, dbexpressions );
  ADD_DEBUG( Sqrt( kInGeV ), dbexpressions );
  return Sqrt( kInGeV );
}

Expression AmpGen::BlattWeisskopf_Norm( const Expression& z2, const Expression& z02, unsigned int L )
{
  if ( L == 0 )
    return Expression( Constant( 1 ) );
  else if ( L == 1 )
    return ( 1 + z02 ) / ( 1 + z2 );
  else if ( L == 2 )
    return ( z02 * z02 + 3 * z02 + 9 ) / ( z2 * z2 + 3 * z2 + 9 );
  else {
    ERROR( "Cannot understand L > 2" );
    return Expression( Constant( 0 ) );
  }
}

Expression AmpGen::BlattWeisskopf( const Expression& z2, unsigned int L )
{
  if ( L == 0 )
    return Expression( Constant( 1 ) );
  else if ( L == 1 )
    return 2 * z2 / ( 1 + z2 );
  else if ( L == 2 )
    return ( 13 * z2 * z2 ) / ( z2 * z2 + 3 * z2 + 9 );
  else {
    ERROR( "Cannot understand L > 2" );
    return Expression( Constant( 0 ) );
  }
}

Expression AmpGen::BL( const Expression& s, const Expression& s0, const Expression& s1, const Expression& s2,
                       const Expression& radius, const unsigned int& L )
{
  const Expression q2  = Abs( Q2( s, s1, s2 ) );
  const Expression q20 = Abs( Q2( s0, s1, s2 ) );
  return Sqrt( BlattWeisskopf( q2 * radius * radius, L ) / BlattWeisskopf( q20 * radius * radius, L ) );
}

std::vector<Expression> AmpGen::parameterVector( const std::string& name, const unsigned int& nParam )
{
  std::vector<Expression> returnVector;
  for ( unsigned int i = 0; i < nParam; ++i ) returnVector.push_back( Parameter( name + std::to_string( i ) ) );
  return returnVector;
}
Expression AmpGen::width( const Expression& s, const Expression& s1, const Expression& s2, const Expression& mass,
                          const Expression& width, const Expression& radius, unsigned int L,
                          DebugSymbols* dbexpressions )
{
  auto q2v = make_cse( Q2(s,s1,s2) );
  const Expression q2  = Ternary( q2v > 0, q2v, 0 );
  const Expression q20 = Abs( Q2( mass * mass, s1, s2 ) );
  Expression BF        = BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, L );
  Expression qr        = Sqrt( q2 / q20 );
  for ( unsigned int i = 0; i < L; ++i ) qr = qr * ( q2 / q20 );

  const Expression mreco = Sqrt( s );
  const Expression mr    = mass / mreco;

  ADD_DEBUG( q2, dbexpressions );
  ADD_DEBUG( q20, dbexpressions );
  ADD_DEBUG( Sqrt( q2 / q20 ), dbexpressions );
  ADD_DEBUG( BF, dbexpressions );
  ADD_DEBUG( qr, dbexpressions );
  ADD_DEBUG( mr, dbexpressions );
  ADD_DEBUG( qr * mr, dbexpressions );
  ADD_DEBUG( Sqrt( q2 ) / mreco, dbexpressions );
  ADD_DEBUG( Sqrt( q20 ) / mass, dbexpressions );
  ADD_DEBUG( width * BF * qr * mr, dbexpressions );

  return width * BF * qr * mr;
}

Expression kFactor( const Expression& mass, const Expression& width, DebugSymbols* dbexpressions )
{
  const Expression massInGeV  = mass / GeV;
  const Expression widthInGeV = width / GeV;
  const Expression gammaInGeV = massInGeV * Sqrt( massInGeV * massInGeV + widthInGeV * widthInGeV );
  const Expression kInGeV =
      2 * sqrt( 2 ) * massInGeV * widthInGeV * gammaInGeV / ( M_PI * Sqrt( massInGeV * massInGeV + gammaInGeV ) );
  ADD_DEBUG( gammaInGeV, dbexpressions );
  ADD_DEBUG( Sqrt( kInGeV ), dbexpressions );
  return Sqrt( kInGeV );
}

//// Lineshape factory declarations ////

template <>
Factory<AmpGen::ILineshape>* Factory<AmpGen::ILineshape>::gImpl = nullptr;

bool LineshapeFactory::isLineshape( const std::string& lineshape )
{

  size_t pos = lineshape.find( "." );
  if ( pos == std::string::npos )
    return get( lineshape, true ) != nullptr;
  else
    return get( lineshape.substr( 0, pos ), true ) != nullptr;
}

Expression LineshapeFactory::getLineshape( const std::string& lineshape, const Expression& s, const Expression& s1,
                                           const Expression& s2, const std::string& particleName, const unsigned int& L,
                                           DebugSymbols* dbexpressions )
{
  size_t pos = lineshape.find( "." );

  if ( pos == std::string::npos ) {
    auto it = get( lineshape );
    if ( !it ) ERROR( "Lineshape : " << lineshape << " not found" );
    return it->get( s, s1, s2, particleName, L, "", dbexpressions );
  } else {
    return get( lineshape.substr( 0, pos ) )
        ->get( s, s1, s2, particleName, L, lineshape.substr( pos + 1 ), dbexpressions );
  }
}

Expression LineshapeFactory::getGenericShape( const std::string& lineshape, const std::vector<Tensor>& p,
                                              const std::string& particleName, const unsigned int& L,
                                              DebugSymbols* dbexpressions )
{
  size_t pos = lineshape.find( "." );

  if ( pos == std::string::npos ) {
    auto it = get( lineshape );
    if ( !it ) ERROR( "Lineshape : " << lineshape << " not found" );
    return it->get( p, particleName, L, "", dbexpressions );
  } else {
    return get( lineshape.substr( 0, pos ) )->get( p, particleName, L, lineshape.substr( pos + 1 ), dbexpressions );
  }
}

Expression AmpGen::pol( const Expression& X, const std::vector<Expression>& p )
{
  Expression F = 0;
  Expression L = 1;
  for ( auto& ip : p ) {
    F = F + ip * L;
    L = L * X;
  }
  return F;
}
