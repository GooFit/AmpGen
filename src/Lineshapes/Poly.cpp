#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

DEFINE_LINESHAPE( Poly )
{
  unsigned int degree = NamedParameter<unsigned int>( lineshapeModifier + "::Degree" );
  auto params = parameterVector( lineshapeModifier + "_c", degree );
  return pol(s, params);
}
