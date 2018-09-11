#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"

using namespace AmpGen;

DEFINE_LINESHAPE( Gaussian )
{
  Expression mu           = Parameter( lineshapeModifier + "_mean" );
  Expression sigma        = Parameter( lineshapeModifier + "_sigma" );
  const Expression sInGeV = s / ( 1000 * 1000 );
  const Expression d      = sInGeV - mu;

  return Exp( -d * d / ( 2 * sigma * sigma ) ) ;
}
