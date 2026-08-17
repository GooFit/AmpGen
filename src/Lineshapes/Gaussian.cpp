#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"

using namespace AmpGen;

DEFINE_LINESHAPE(Gaussian) {
  Expression mu = Variable(lineshapeModifier + "_mean");
  Expression sigma = Variable(lineshapeModifier + "_sigma");
  const Expression d = s - mu;
  return fcn::exp(-d * d / (2 * sigma * sigma));
}
