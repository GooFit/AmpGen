#ifndef AMPGEN_WIGNER_H
#define AMPGEN_WIGNER_H
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Transform.h"

namespace AmpGen {
  class Particle;

  Expression wigner_d( const Expression& cb, const double& j, const double& m, const double& n );
  Expression wigner_D( const std::pair<Expression, Expression>& P, const double& J, const double& lA, const double& lB, DebugSymbols* db);
  /** @ingroup Vertices function CG  
    Calculates the Clebsch-Gordan coefficient for (j1 m1 j2 m2 | J M), the expansion
    coefficients in  
    */
  double CG( const double& j1,
      const double& m1,
      const double& j2, 
      const double& m2,
      const double& J,
      const double& M );

  /** @ingroup Vertices function wickTransform 
    Generates a wick transform sequence that aligns tensor P (four-vector) to the +/- ve z-axis, then boosts to the rest frame. 
    The mass may be seperately specified. The parameter ve specifies whether the initial Euler rotation is to the +/- z-axis. 
    In the case where ve =-1, a second rotation is applied about the x-axis that aligns P to the +ve z-axis. 
    This ensures that singly and doubly primed helicity frames remain orthonormal.
    */
  TransformSequence wickTransform(const Tensor& P, const Expression& M, const int& ve =1, const bool& handleZeroCase = false);   

  Expression helicityAmplitude(const Particle& particle, TransformSequence& parentFrame, const double& Mz,    DebugSymbols* db , const int sgn=1, std::map<const Particle*, TransformSequence>* cacheptr = nullptr); 
  Tensor basisSpinor(const int& polState, const int& id);
  Tensor basisVector(const int& polState);

  struct LS {
    double factor;
    double cg1;
    double cg2;
    double p; 
    double m1;
    double m2;
  };

  std::vector<LS> calculate_recoupling_constants( 
      const double& J, 
      const double& M,
      const double& L, 
      const double& S,
      const double& j1,
      const double& j2 );
}

#endif
