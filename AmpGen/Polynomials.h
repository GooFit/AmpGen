#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"


#ifndef POLYNOMIALS_H_
#define POLYNOMIALS_H_



namespace AmpGen { 

Expression CPSinPoly(Expression x, Expression y, unsigned int mu, unsigned int lambda, int CP){
  auto M00 = fcn::sin(M_PI * lambda * x);
  auto M01 = fcn::sin(M_PI * lambda * y);
  auto M10 = fcn::sin(M_PI * mu * x);
  auto M11 = fcn::sin(M_PI * mu * y);
  return M00 * M11 + CP*M10*M01;
}

Expression chebychev(Expression x, unsigned int order){
    Expression output = 0;
    if (order == 0){
        output = 1;
    }
    else if (order == 1){
        output = x;
    }
    else if (order == 2){
        output = 2 * x * x - 1;
    }
    else {
        output = 2 * x * chebychev(x, order -1) - chebychev(x, order - 2);
    }
    return output;
}


Expression legendre(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x;
    }
    else {
        output = (2 * order - 1)/order * x * legendre(x, order - 1) - (order - 1)/order * legendre(x, order - 2);
    }
    return output;
}
    
Expression bessel(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x + 1;
    }
    else{
        output = (2 * order - 1) * x * bessel(x, order - 1) + bessel(x, order - 2);
        
    }
    return output;
}

Expression laguerre(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output= 1 - x;
    }
    else{
        output = ( (2 * order  - 1 - x) * laguerre(x, order - 1) - (order - 1) * laguerre(x, order - 2))/order;
    }
    return output;
}
}