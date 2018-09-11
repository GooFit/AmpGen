#ifndef AMPGEN_DIRACMATRICES_H
#define AMPGEN_DIRACMATRICES_H

#include "AmpGen/Tensor.h"

namespace AmpGen
{
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
