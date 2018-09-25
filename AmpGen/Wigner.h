#ifndef AMPGEN_WIGNER_H
#define AMPGEN_WIGNER_H
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"

namespace AmpGen {
  Expression wigner_d( const Expression& beta, const double& j, const double& m, const double& n );
}

#endif
