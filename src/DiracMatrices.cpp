#include "AmpGen/DiracMatrices.h"

#include <stddef.h>
#include <initializer_list>

#include "AmpGen/Tensor.h"

using namespace std::complex_literals;
using namespace AmpGen;

extern const AmpGen::Expression  AmpGen::I = Constant(0, 1);
extern const AmpGen::Expression  AmpGen::Z = Constant(0);

extern const std::array<AmpGen::Tensor,5>     AmpGen::Gamma  ( { 
    Tensor ({ 0, 0, 0, 1, 
              0, 0, 1, 0, 
              0,-1, 0, 0, 
             -1, 0, 0, 0} , {4, 4}),
    Tensor ({ Z, Z, Z,-I, 
              Z, Z, I, Z, 
              Z, I, Z, Z, 
             -I, Z, Z, Z} , {4, 4}),
    Tensor ({ 0, 0, 1, 0, 
              0, 0, 0,-1, 
             -1, 0, 0, 0,  
              0, 1, 0, 0} , {4, 4}),
    Tensor ({ 1, 0, 0, 0, 
              0, 1, 0, 0, 
              0, 0,-1, 0, 
              0, 0, 0,-1} , {4, 4}),
    Tensor ({ 0, 0, 1, 0, 
              0, 0, 0, 1, 
              1, 0, 0, 0, 
              0, 1, 0, 0} , {4, 4})} );

extern const std::array<AmpGen::Tensor,3> AmpGen::Sigma ( 
    {
      Tensor({0,1,1,0},{2,2}),
      Tensor({Z,-I,Z,I},{2,2}),
      Tensor({1,0,0,-1},{2,2})
    } );

extern const std::array<AmpGen::Tensor,3>     AmpGen::S03  ( { 
    Tensor ({ 0, 0, 0, 0,  
              0, 0,-1, 0,  
              0, 1, 0, 0,
              0, 0, 0, 0 } , Tensor::dim(4,4) ),
    Tensor ({ 0, 0, 1, 0,  
              0, 0, 0, 0,  
             -1, 0, 0, 0,
              0, 0, 0, 0 } , Tensor::dim(4,4) ),
    Tensor ({ 0,-1, 0, 0,  
              1, 0, 0, 0,  
              0, 0, 0, 0,
              0, 0, 0, 0 } , Tensor::dim(4,4)  ) } );
