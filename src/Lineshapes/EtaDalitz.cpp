#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/Tensor.h"

using namespace AmpGen;

DEFINE_GENERIC_SHAPE( EtaDalitz )
{

  Tensor P( std::vector<size_t>( {4} ) );
  for ( auto& ip : p ) P = P + ip;

  Expression T0 = Sqrt( dot( p[0], P ) * dot( p[0], P ) / dot( P, P ) ) - Sqrt( dot( p[0], p[0] ) );
  Expression T1 = Sqrt( dot( p[1], P ) * dot( p[1], P ) / dot( P, P ) ) - Sqrt( dot( p[1], p[1] ) );
  Expression T2 = Sqrt( dot( p[2], P ) * dot( p[2], P ) / dot( P, P ) ) - Sqrt( dot( p[2], p[2] ) );

  Expression Q   = T0 + T1 + T2;
  Expression y   = 3.0 * T0 / Q - 1.0;
  Expression z   = 1.0 - 1.07 * y;
  Expression amp = Ternary( z > 0.0, Sqrt( z ), 1 ) ; 

  if ( lineshapeModifier != "" ) {
    amp = amp * LineshapeFactory::getLineshape( lineshapeModifier, dot( P, P ), dot( p[0], p[0] ), dot( p[1], p[1] ),
                                                particleName, L, dbexpressions );
  }
  return amp;
}
