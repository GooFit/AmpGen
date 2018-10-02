#ifndef AMPGEN_WIGNER_H
#define AMPGEN_WIGNER_H
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"

namespace AmpGen {
  Expression wigner_d( const Expression& cb, const double& j, const double& m, const double& n );
  Expression wigner_D( const Tensor& P, const double& J, const double& lA, const double& lB, const double& lC ); 

  double CG( const double& j1,
    const double& m1,
    const double& j2, 
    const double& m2,
    const double& J,
    const double& M );

  Tensor rotationMatrix( const Tensor& p );
  Tensor helicityTransformMatrix( const Tensor& P, const Expression& E, const Expression& M, const double& ve =1); 
  
}

#endif
