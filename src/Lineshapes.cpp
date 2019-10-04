#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Units.h"
#include <stddef.h>
#include <cmath>
#include <ostream>

using namespace AmpGen;
using namespace AmpGen::fcn;

Expression AmpGen::Q2( const Expression& s, const Expression& M1sq, const Expression& M2sq )
{
  const Expression num = s -  2*M1sq - 2*M2sq  + ( M1sq - M2sq ) * ( M1sq - M2sq ) / s ;
  return num/4.;
}

Expression AmpGen::kFactor( const Expression& mass, const Expression& width, DebugSymbols* dbexpressions )
{
  const Expression gamma = mass * sqrt( mass * mass + width * width );
  const Expression k = (2 * sqrt(2) / M_PI) * mass * width * gamma / sqrt( mass*mass + gamma );
  ADD_DEBUG( gamma  , dbexpressions );
  ADD_DEBUG( sqrt(k), dbexpressions );
  return sqrt(k);
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

Expression AmpGen::BL( const Expression& s, 
                       const Expression& s0, 
                       const Expression& s1, 
                       const Expression& s2,
                       const Expression& radius, 
                       const unsigned int& L )
{
  const Expression q2  = abs( Q2( s, s1, s2 ) );
  const Expression q20 = abs( Q2( s0, s1, s2 ) );
  return sqrt( BlattWeisskopf( q2 * radius * radius, L ) / BlattWeisskopf( q20 * radius * radius, L ) );
}

std::vector<Expression> AmpGen::parameterVector( const std::string& name, const size_t& nParam )
{
  std::vector<Expression> returnVector;
  for (size_t i = 0; i != nParam; ++i ) 
    returnVector.emplace_back(Parameter(name + std::to_string(i)));
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

bool Lineshape::Factory::isLineshape( const std::string& lineshape )
{

  size_t pos = lineshape.find( "." );
  if ( pos == std::string::npos )
    return AmpGen::Factory<Lineshape::Base>::get(lineshape, true) != nullptr;
  else
    return AmpGen::Factory<Lineshape::Base>::get(lineshape.substr(0, pos), true ) != nullptr;
}

Expression Lineshape::Factory::get( const std::string& lineshape, const Expression& s, const Expression& s1,
                                           const Expression& s2, const std::string& particleName, const unsigned int& L,
                                           DebugSymbols* dbexpressions )
{
  size_t pos = lineshape.find( "." );

  if ( pos == std::string::npos ) {
    auto it = AmpGen::Factory<Lineshape::Base>::get( lineshape );
    if ( !it ) ERROR( "Lineshape : " << lineshape << " not found" );
    return it->get( s, s1, s2, particleName, L, "", dbexpressions );
  } else {
    return AmpGen::Factory<Lineshape::Base>::get(lineshape.substr(0, pos))->get( s, s1, s2, particleName, L, lineshape.substr( pos + 1 ), dbexpressions );
  }
}

Expression Lineshape::Factory::get(const std::string& lineshape, const AmpGen::Particle& p,
                                   DebugSymbols* dbexpressions )
{
  size_t pos = lineshape.find( "." );
  
  if ( pos == std::string::npos ) {
    auto it = AmpGen::Factory<Lineshape::Base>::get( lineshape );
    if ( !it ) ERROR( "Lineshape : " << lineshape << " not found" );
    return it->get(p, "", dbexpressions );
  } else {
    return AmpGen::Factory<Lineshape::Base>::get( lineshape.substr( 0, pos ) )->get(p, lineshape.substr( pos + 1 ), dbexpressions );
  }
}

Expression AmpGen::pol( const Expression& X, const std::vector<Expression>& p )
{
  Expression F = 0;
  Expression L = 1;
  for ( auto& ip : p ) {
    F = F + ip * L;
    L = L * X;
/*     std::cout<<"F = "<<F<<"\n";
    std::cout<<"ip = "<<ip<<"\n";
    std::cout<<"X = "<<X<<"\n";
    std::cout<<"L = "<<L<<"\n";
    */
  }
  return F;
}
