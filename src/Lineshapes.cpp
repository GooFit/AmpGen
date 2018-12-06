#include "AmpGen/Lineshapes.h"

#include <stddef.h>
#include <cmath>
#include <ostream>

#include "AmpGen/MsgService.h"
#include "AmpGen/Units.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

template <>
Factory<AmpGen::ILineshape>* Factory<AmpGen::ILineshape>::gImpl = nullptr;

Expression AmpGen::Q2( const Expression& s, const Expression& M1sq, const Expression& M2sq )
{
  const Expression num = s -  2*M1sq - 2*M2sq  + ( M1sq - M2sq ) * ( M1sq - M2sq ) / s ;
  return num/4.;
}

Expression AmpGen::kFactor( const Expression& mass, const Expression& width, DebugSymbols* dbexpressions )
{
  const Expression massInGeV  = mass / GeV;
  const Expression widthInGeV = width / GeV;
  const Expression gammaInGeV = massInGeV * sqrt( massInGeV * massInGeV + widthInGeV * widthInGeV );
  const Expression kInGeV =
      2 * sqrt( 2 ) * massInGeV * widthInGeV * gammaInGeV / ( M_PI * sqrt( massInGeV * massInGeV + gammaInGeV ) );
  ADD_DEBUG( gammaInGeV  , dbexpressions );
  ADD_DEBUG( sqrt(kInGeV), dbexpressions );
  return sqrt( kInGeV );
}

Expression AmpGen::BlattWeisskopf_Norm( const Expression& z, const Expression& z0, unsigned int L )
{
  switch( L ) {
    case 0: return 1;
    case 1: return (1+z0) / (1+z);
    case 2: return (z0*z0 + 3*z0 + 9 ) / (z*z + 3*z + 9);
    case 3: return (z0*fpow(z0-15,2) + 9*fpow( 2*z0 -5,2 ) ) / (z*fpow(z-15,2) + 9*fpow( 2*z -5,2 ) );
    case 4: return (fpow(z0*z0-45*z0+105, 2)+25*z*fpow(2*z-21,2)) / (fpow(z*z-45*z+105, 2)+25*z*fpow(2*z-21,2) );
    default : 
      ERROR("Barrier factors not for implemented for L>4");
      return 0;
  }
}

Expression AmpGen::BlattWeisskopf( const Expression& z, unsigned int L )
{
  switch( L ) {
    case 0: return 1;
    case 1: return 2    *fpow(z,1) / (1+z);
    case 2: return 13   *fpow(z,2) / (z*z + 3*z + 9);
    case 3: return 277  *fpow(z,3) / (z*fpow(z-15,2) + 9*fpow( 2*z -5,2 ) );
    case 4: return 12746*fpow(z,4) / (fpow(z*z - 45*z + 105, 2) + 25*z*fpow(2*z-21,2) );
    default : 
      ERROR("Barrier factors not for implemented for L>4");
      return 0;
  }
}

Expression AmpGen::BL( const Expression& s, const Expression& s0, const Expression& s1, const Expression& s2,
                       const Expression& radius, const unsigned int& L )
{
  const Expression q2  = abs( Q2( s, s1, s2 ) );
  const Expression q20 = abs( Q2( s0, s1, s2 ) );
  return sqrt( BlattWeisskopf( q2 * radius * radius, L ) / BlattWeisskopf( q20 * radius * radius, L ) );
}

std::vector<Expression> AmpGen::parameterVector( const std::string& name, const unsigned int& nParam )
{
  std::vector<Expression> returnVector;
  for ( unsigned int i = 0; i < nParam; ++i ) 
    returnVector.push_back( Parameter( name + std::to_string( i ) ) );
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
  Expression qr        = sqrt( q2 / q20 ) * fpow( q2/q20, L );

  const Expression mreco = isqrt( s );
  const Expression mr    = mass * mreco;

  ADD_DEBUG(q2            , dbexpressions);
  ADD_DEBUG(q20           , dbexpressions);
  ADD_DEBUG(sqrt(q2/q20)  , dbexpressions);
  ADD_DEBUG(BF            , dbexpressions);
  ADD_DEBUG(qr            , dbexpressions);
  ADD_DEBUG(mr            , dbexpressions);
  ADD_DEBUG(qr*mr         , dbexpressions);
  ADD_DEBUG(sqrt(q2)/mreco, dbexpressions);
  ADD_DEBUG(sqrt(q20)/mass, dbexpressions);
  ADD_DEBUG(width*BF*qr*mr, dbexpressions);

  return width * BF * qr * mr;
}

bool Lineshape::LineshapeFactory::isLineshape( const std::string& lineshape )
{

  size_t pos = lineshape.find( "." );
  if ( pos == std::string::npos )
    return get( lineshape, true ) != nullptr;
  else
    return get( lineshape.substr( 0, pos ), true ) != nullptr;
}

Expression Lineshape::LineshapeFactory::getLineshape( const std::string& lineshape, const Expression& s, const Expression& s1,
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

Expression Lineshape::LineshapeFactory::getGenericShape( const std::string& lineshape, const std::vector<Tensor>& p,
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
