#include "AmpGen/Expression.h"
#include "AmpGen/Particle.h"

namespace AmpGen { 

  namespace Lineshape {
    /** @ingroup Lineshapes class CoupledChannel 
      @brief Description of a resonance that decays to multiple two and three-body final states. */
    DECLARE_LINESHAPE( CoupledChannel );
  }
  Expression phaseSpace(const Expression& s, const Particle& p, const size_t& l);
}
