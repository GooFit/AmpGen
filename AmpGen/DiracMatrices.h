#ifndef AMPGEN_DIRACMATRICES_H
#define AMPGEN_DIRACMATRICES_H

#include <array>

#include "AmpGen/Tensor.h"
#include "AmpGen/Expression.h"

namespace AmpGen
{
  extern const Expression I;
  extern const Expression Z;
  extern const std::array<Tensor,5> Gamma;
  extern const std::array<Tensor,3> Sigma;
  extern const std::array<Tensor,3> S03;
  extern const std::array<Tensor,8> SU3;
} // namespace AmpGen

#endif
