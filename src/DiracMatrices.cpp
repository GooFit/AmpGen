#include "AmpGen/DiracMatrices.h"

#include <stddef.h>
#include <initializer_list>

#include "AmpGen/Tensor.h"

using namespace std::complex_literals;
using namespace AmpGen;

extern const AmpGen::Expression  AmpGen::I = Constant(0, 1);
extern const AmpGen::Expression  AmpGen::Z = Constant(0);

extern const std::array<AmpGen::Tensor,5>     AmpGen::Gamma( { 
    Tensor({ 0, 0, 0, 1, 0, 0, 1, 0, 0,-1, 0, 0,-1, 0, 0, 0}, Tensor::dim(4,4)),
    Tensor({ Z, Z, Z,-I, Z, Z, I, Z, Z, I, Z, Z,-I, Z, Z, Z}, Tensor::dim(4,4)),
    Tensor({ 0, 0, 1, 0, 0, 0, 0,-1,-1, 0, 0, 0, 0, 1, 0, 0}, Tensor::dim(4,4)),
    Tensor({ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1}, Tensor::dim(4,4)),
    Tensor({ 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0}, Tensor::dim(4,4))} );

extern const std::array<AmpGen::Tensor,3> AmpGen::Sigma( {
    Tensor({ 0, 1, 1, 0}, Tensor::dim(2,2)),
    Tensor({ Z,-I, I, Z}, Tensor::dim(2,2)),
    Tensor({ 1, 0, 0,-1}, Tensor::dim(2,2))} );

extern const std::array<AmpGen::Tensor,3>     AmpGen::S03  ( { 
    Tensor({ 0, 0, 0, 0, 0, 0,-1, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, Tensor::dim(4,4) ),
    Tensor({ 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0 }, Tensor::dim(4,4) ),
    Tensor({ 0,-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, Tensor::dim(4,4)  ) } );

extern const std::array<AmpGen::Tensor,8>     AmpGen::SU3 ( {
    Tensor({ 0, 1, 0, 1, 0, 0, 0, 0, 0}, Tensor::dim(3,3)),
    Tensor({ Z,-I, Z, I, Z, Z, Z, Z, Z}, Tensor::dim(3,3)),
    Tensor({ 1, 0, 0, 0,-1, 0, 0, 0, 0}, Tensor::dim(3,3)),
    Tensor({ 0, 0, 1, 0, 0, 0, 1, 0, 0}, Tensor::dim(3,3)),
    Tensor({ Z, Z,-I, Z, Z, Z, I, Z, Z}, Tensor::dim(3,3)),
    Tensor({ 0, 0, 0, 0, 0, 1, 0, 1, 0}, Tensor::dim(3,3)),
    Tensor({ Z, Z, Z, Z, Z,-I, Z, I, Z}, Tensor::dim(3,3)),
    Tensor({ 1./sqrt(3), 0., 0., 0.,1./sqrt(3.), 0., 0., 0., -2./sqrt(3.)}, Tensor::dim(3,3)) });
