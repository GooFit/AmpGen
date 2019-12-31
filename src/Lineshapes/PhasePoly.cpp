#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/NamedParameter.h"

#include <cmath>

using namespace AmpGen;
using namespace AmpGen::fcn;

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

DEFINE_GENERIC_SHAPE( PhasePoly )
{
    unsigned int degree = NamedParameter<unsigned int>( lineshapeModifier + "::Degree" ) ;
    //vector of vectors 
    std::vector<std::vector<Expression> > param;
    unsigned int i=0;
    bool debug = NamedParameter<bool>("PhasePoly::Debug", false);
    std::string polyType= NamedParameter<std::string>("PhasePoly::PolyType", "simple");

    if (debug){
            INFO("2D Polynomial of the form P(x,y), used to correct the fitted phase of D->KsPiPi decays ");
            INFO("Using "<<degree<<" order "<< polyType <<" polynomial");
            //Degree for our 2D polynomial - N
           

    }

    for (int i=0; i<degree+1; i++){
        //auto param_i = parameterVector(lineshapeModifier + "::C" + std::to_string(i), degree - i);
        std::vector<Expression> param_i;
        for (int j=0; j < degree+1-i; j++){
            auto param_ij = Parameter(lineshapeModifier+"::C"+std::to_string(i)+std::to_string(j));
            param_i.push_back(param_ij);
        }
        if (debug){
            INFO(lineshapeModifier+"::C" + std::to_string(i) + ", " + std::to_string(degree-i));
            int j=0;
            for (auto paramij : param_i){
                INFO("C"<<i<<j<<" = "<<paramij);
                j++;
            }
        }
       
        param.push_back(param_i);
    }
    

    //4D momentum tensor for 

    //Tensor P (Tensor::dim(4));
    //for ( auto& ip : p ) P = P + ip;
    //For simplicity we will use x and y as the input for the polynomial not m^2_+, m^2_-!
  //  Expression x = dot(p[0] + p[1], p[0] + p[1]);
//    Expression y = dot(p[0] + p[2], p[0] + p[2]);
//    Expression x = dot(p[0], p[0]);
//    Expression y = x;
    auto pp = *p.daughter("pi+");
    auto pm = *p.daughter("pi-");
    auto ks = *p.daughter("K0S0");
    auto pD = ks.P() + pp.P() + pm.P();
    auto mD = sqrt(dot(pD,pD));
    auto mp = sqrt(dot(pp.P(), pp.P()));
    auto mm = sqrt(dot(pm.P(), pm.P()));
    auto mK = sqrt(dot(ks.P(), ks.P()));


    //Tensor P (Tensor::dim(4));
    //for ( auto& ip : p ) P = P + ip;
    //For simplicity we will use x and y as the input for the polynomial not m^2_+, m^2_-!

    Expression xmin = pow(mp + mK, 2);
    Expression xmax = pow(mD - mm, 2);
    Expression x0 = (xmax + xmin)/2;
    Expression ymin = pow(mp + mK, 2);
    Expression ymax = pow(mD - mp, 2);
    Expression y0 = (ymax + ymin)/2;
    Expression x = dot(ks.P() + pp.P(), ks.P() + pp.P())/x0 - 1;
    Expression y = dot(ks.P() + pm.P(), ks.P() + pm.P())/y0 - 1;
    Expression sum =0;
    //For a 2D polynomial we take the x projections so V_i = (c_i0, ci1,.., c_im), where m = N+1 - i
    //Then we do the y polynomial for each x^i, then sum all of the y polynomials.
    //i.e. sum_0 = c_00 + c_01 y + c_02 y^2 ...
    //     sum_1 = c_10 x + c_11 xy + c12 x y^2 + ...
    //     sum_2 = c_20 x^2 + c_21 x^2 y + c_22 x^2 y^2 + ...

    i=0;
    while (i!=degree){
        Expression sum_i=0;
        for (unsigned int j=0; j<param[i].size(); j++){
            if (debug){
            INFO("x = "<<x());
            INFO("y = "<<y());
            INFO("param["<<i<<"]["<<j<<"] ="<<param[i][j]);
            INFO("Type = "<<polyType) ;
            }

            if (polyType=="simple"){
                sum_i = sum_i +  param[i][j] * pow(x, i) * pow(y, j);
            }
            else if (polyType=="chebychev"){
               sum_i = sum_i + param[i][j] * chebychev(x, i) * chebychev(y, j); 
            }
            else if (polyType=="legendre"){
                sum_i = sum_i + param[i][j] * legendre(x, i) * legendre(y, j);
            }
            else if (polyType=="laguerre"){
                sum_i = sum_i + param[i][j] * laguerre(x, i) * laguerre(y, j);
            }
            if (debug){
            INFO("sum_"<<i<<" = "<<sum_i());
                    }
            
        }
        i++;
        sum = sum + sum_i;
    }
   // Expression amp = exp(std::complex<double>(0,1) * sum);
    //Expression amp = sum;
    auto amp = exp(Constant(0,1) * sum);
   // std::cout<<"amp = "<<amp<<"\n";
    return amp;
}


