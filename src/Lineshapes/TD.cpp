#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/ParticleProperties.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

DEFINE_GENERIC_SHAPE( TD )
{
  Expression tau = Parameter(p.name() +"_decayTime");   
  ADD_DEBUG( tau, dbexpressions );
  ADD_DEBUG( p.props()->lifetime(), dbexpressions );
  return tau / ( 2 * p.props()->lifetime() ); 
}
