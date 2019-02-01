#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/Tensor.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

DEFINE_GENERIC_SHAPE( EtaDalitz )
{
  WARNING("Empirical expression; eta0[EtaDalitz.NonRelBW]{pi+,pi-,pi0} with pi0 in the 3rd order.");
  Tensor P( Tensor::dim(4) );
  for ( auto& ip : p ) P = P + ip;

  Expression T0 = sqrt( dot( p[0], P ) * dot( p[0], P ) / dot( P, P ) ) - sqrt( dot( p[0], p[0] ) );
  Expression T1 = sqrt( dot( p[1], P ) * dot( p[1], P ) / dot( P, P ) ) - sqrt( dot( p[1], p[1] ) );
  Expression T2 = sqrt( dot( p[2], P ) * dot( p[2], P ) / dot( P, P ) ) - sqrt( dot( p[2], p[2] ) );

  Expression Q   = T0 + T1 + T2;
  Expression y   = 3.0 * T2 / Q - 1.0;
  Expression z   = 1.0 - 1.07 * y;
  Expression amp = Ternary( z > 0.0, sqrt(z), 1 ); 

  if ( lineshapeModifier != "" ){
    amp = amp * LineshapeFactory::getGenericShape(lineshapeModifier, s, p, particleName, L, dbexpressions);
  }
  return amp;
}

