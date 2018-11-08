#ifndef AMPGEN_DIRACMATRICES_H
#define AMPGEN_DIRACMATRICES_H

#include <array>

#include "AmpGen/Tensor.h"
#include "AmpGen/Expression.h"

namespace AmpGen
{
class Tensor;

  extern const Expression I;
  extern const Expression Z;
  namespace Dirac
  {
    extern const std::array<Tensor, 5> Gamma;
  }
  namespace Weyl 
  {
    extern const std::array<Tensor, 5> Gamma;
  }
  extern const std::array<Tensor,3> Sigma;

} // namespace AmpGen

#endif
