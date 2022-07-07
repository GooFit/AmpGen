#ifndef AMPGEN_DIRACMATRICES_H
#define AMPGEN_DIRACMATRICES_H

#include <array>

#include "AmpGen/Tensor.h"
#include "AmpGen/Expression.h"

namespace AmpGen
{
  extern const Expression I;                // Imaginary unit
  extern const Expression Z;                // Zero
  extern const std::array<Tensor,5> Gamma;  // Gamma matrices, orderered weirdly due to ROOT convention, so Gamma[{0,1,2}] = gamma[{1,2,3}] and Gamma[3] = gamma[0] 
                                            // where gamma are the gamma matrices in the more conventional ordering. 
  extern const std::array<Tensor,3> Sigma;  // The pauli matrices 
  extern const std::array<Tensor,3> Sigma4; // The pauli matrices acting on bispinors  
  extern const std::array<Tensor,3> SO3;    // basis elements for the vector reprsentation of SO(3)
  extern const std::array<Tensor,8> SU3;    // basis elements for the vector representation of SU(3), i.e. the Gell-Mann matrices
} // namespace AmpGen

#endif
