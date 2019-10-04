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
  auto pp = *p.daughter("pi+");
  auto pm = *p.daughter("pi-");
  auto p0 = *p.daughter("pi0");
  auto  P = p.P();

  Expression T0 = dot(pp.P(), P) - p.mass() * pp.mass();
  Expression T1 = dot(pm.P(), P) - p.mass() * pm.mass();
  Expression T2 = dot(p0.P(), P) - p.mass() * p0.mass();

  Expression Q   = T0 + T1 + T2;
  Expression y   = 3.0 * T2 / Q - 1.0;
  Expression z   = 1.0 - 1.07 * y;
  Expression amp = Ternary( z > 0.0, sqrt(z), 1 ); 

  if ( lineshapeModifier != "" ){
    amp = amp * Lineshape::Factory::get(lineshapeModifier, p, dbexpressions);
  }
  return amp;
}

