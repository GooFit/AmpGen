
#ifndef POLYNOMIALS_H_
#define POLYNOMIALS_H_


#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"



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

    if (order==0){
        return 1;
    }
    else if (order==1){
        return x;
    }
    else {
        return (2 * order - 1)/order * x * legendre(x, order - 1) - (order - 1)/order * legendre(x, order - 2);
    }

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




Expression customPoly(Expression x, Expression y, unsigned int i, unsigned int j){
    auto c00 = legendre(x, 0) * legendre(y, 1);
    auto c01 = legendre(x, 0) * legendre(y, 3);
    auto c10 = legendre(x, 1) * legendre(y, 1);
    auto c11 = legendre(x, 1) * legendre(y, 3);
    auto c02 = legendre(x, 0) * legendre(y, 5);
    auto c20 = legendre(x, 2) * legendre(y, 1);

    //real_t M [3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
//    real_t M [3][3] = {{1,-1,1}, {-1,1,-1}, {1,-1,1}};
   //   real_t M [3][3] = {{-1,1,1}, {0,1,1}, {1,0,1}};
//    real_t M [3][3] = {{-1,1,1}, {0,1,-1}, {1,0,1}};
    //real_t M [3][3] = {{-1/3,1/3,2/3}, {1/3,2/3,1/3}, {1/3,-1/3,1/3}};
    //real_t M [3][3] = {{1/3,-1/3,1/3}, {1/3,2/3,1/3}, {-1/3,1/3,2/3}};

//delta_3 = 1/8
//S
//    real_t M [3][3] = {{-1, 1, 1 }, {0, (1/16) * (7 + pow(561, 0.5)), (1/16)*(7 - pow(561, 0.5))}, {1, 1, 1}};
//S^-1
//    real_t M [3][3] = {{-1/2, 0, 1/2}, {1/4 - 7/(4 * pow(561, 0.5)), 8/pow(561, 0.5), 1/4 - 7/(4 * pow(561, 0.5)) }, {1/4 + 7/(4 * pow(561, 0.5)), -8/pow(561, 0.5), 1/4 + 7/(4 * pow(561, 0.5))}  };
   
//    real_t M [3][3] = {{1,0,1}, {0,1,-1}, {-1,1,1}};
//    real_t M [3][3] = {{1/3,-1/3,1/3}, {1/3,2/3,1/3}, {-1/3,1/3,2/3}};


//delta_3 = 1/4
//    real_t M[3][3] = {{-1, 1, 1}, {0, (1/8) * (3 + pow(137, 0.5)), (1/8) * (3 - pow(137, 0.5))}, {1,1,1}};
//S
//S^-1
    real_t M[3][3] = {{-1/2, 0, -1/2}, { 1/4 - 3/(4 * pow(137, 0.5)), 4/pow(137, 0.5), 1/4 - 3/(4 * pow(137, 0.5))  }, { 1/4 + 3/(4 * pow(137, 0.5)), -4/pow(137, 0.5), 1/4 + 3/(4 * pow(137, 0.5)) }};
//aS^-1
//    real_t M[3][3] = {{ 1/4 + 3/(4 * pow(137, 0.5)), -4/pow(137, 0.5), 1/4 + 3/(4 * pow(137, 0.5)) } , { 1/4 - 3/(4 * pow(137, 0.5)), 4/pow(137, 0.5), 1/4 - 3/(4 * pow(137, 0.5))  }, {-1/2, 0, 1/2}  };


//Wrong S

//    real_t M [3][3] = {{-1,1,-1}, {0,1,-1}, {1,0,1}};


    if (i==0 && j ==0){
        //return 1/3 * (-c00 + c01 + 2* c10);
        //return c00;
        return M[0][0] * c00 + M[0][1] * c01 + M[0][2] * c10;
    }
    else if (i==0 && j==1){
        //return  (c01 - c10);
        return M[1][0] * c00 + M[1][1] * c01 + M[1][2] * c10;
    }
    else if (i==1 && j==0){
        //return 1/3 * (c00 - c01 + c10);
        return M[2][0] * c00 + M[2][1] * c01 + M[2][2] * c10;
        //return (c01 + c10);
        //return c10;
    }
    if (i==2 && j ==0){
        //return 1/3 * (-c00 + c01 + 2* c10);
        //return c00;
        return M[0][0] * c02 + M[0][1] * c11 + M[0][2] * c20;
    }
    else if (i==1 && j==1){
        //return  (c01 - c10);
        return M[1][0] * c02 + M[1][1] * c11 + M[1][2] * c20;
    }
    else if (i==0 && j==2){
        //return 1/3 * (c00 - c01 + c10);
        return M[2][0] * c02 + M[2][1] * c11 + M[2][2] * c20;
        //return (c01 + c10);
        //return c10;
    }




    else {
        return 0;
    }
    }


}
#endif