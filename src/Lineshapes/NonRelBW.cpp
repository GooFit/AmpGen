#include <cmath>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

DEFINE_LINESHAPE( NonRelBW )
{
  auto props          = ParticlePropertiesList::get( particleName );
  Expression mass     = Parameter( particleName + "_mass", props->mass() );
  Expression width    = Parameter( particleName + "_width", props->width() );
  const Expression BW = sqrt( width / 2.0 / M_PI ) / ( sqrt(s) - mass  - Constant(0,1) * width / 2.0 );
  ADD_DEBUG( sqrt(s), dbexpressions );
  ADD_DEBUG( width, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return BW;
}
