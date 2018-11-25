#include "AmpGen/DiracMatrices.h"

#include <stddef.h>
#include <initializer_list>

#include "AmpGen/Tensor.h"

using namespace AmpGen;

extern const AmpGen::Expression               AmpGen::I      ( Constant(0., 1.) );
extern const AmpGen::Expression               AmpGen::Z      ( Expression( 0 ) );

extern const std::array<AmpGen::Tensor,5>     AmpGen::Dirac::Gamma  ( { 
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

extern const std::array<AmpGen::Tensor,5>     AmpGen::Weyl::Gamma  ( { 
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
    Tensor ({ 0, 0, 1, 0, 
              0, 0, 0, 1, 
              1, 0, 0, 0, 
              0, 1, 0, 0} , {4, 4}),
    Tensor ({-1, 0, 0, 0, 
              0,-1, 0, 0, 
              0, 0, 1, 0,
              0, 0, 0, 1} , {4, 4} ) }  );

extern const std::array<AmpGen::Tensor,3> AmpGen::Sigma ( 
    {
      Tensor({0,1,1,0},{2,2}),
      Tensor({Z,-I,I,Z},{2,2}),
      Tensor({1,0,0,-1},{2,2})
    } );
