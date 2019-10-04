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

DEFINE_GENERIC_SHAPE( PhasePoly )
{
    INFO("2D exponential Polynomial of the form e^(iP(x,y)), used to correct the fitted phase of D->KsPiPi decays ");
    //Degree for our 2D polynomial - N
    std::cout<<"Starting to make Phase Poly\n";
    unsigned int degree = NamedParameter<unsigned int>( lineshapeModifier + "::Degree" );
    //vector of vectors 
    std::vector<std::vector<Expression> > param;
    unsigned int i=0;
    while (i!= degree+1){
        auto param_i = parameterVector(lineshapeModifier + "_c" + std::to_string(i), degree + 1 - i);
        i++;
        param.push_back(param_i);
    }
    bool debug = false;
    //4D momentum tensor for 
    //auto sizeP = p.size();
    //if (sizeP != 3){
     //   INFO("I expect 3 particles!");
     //   return 1;
   // }
   auto pp = *p.daughter("pi+");
   auto pm = *p.daughter("pi+");
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
    Expression x = dot(ks.P() + pp.P(), ks.P() + pp.P())/x0;
    Expression y = dot(ks.P() + pm.P(), ks.P() + pm.P())/y0;
   

    //Tensor P (Tensor::dim(4));
    //for ( auto& ip : p ) P = P + ip;
    //For simplicity we will use x and y as the input for the polynomial not m^2_+, m^2_-!

   
   

  //      Expression x = dot(p[0], p[0]);
//        Expression y = dot(p[0], p[0]);
//    Expression x = dot(p[0], p[0]);
//    Expression y = x;
   
    Expression sum =0;
    //For a 2D polynomial we take the x projections so V_i = (c_i0, ci1,.., c_im), where m = N+1 - i
    //Then we do the y polynomial for each x^i, then sum all of the y polynomials.
    //i.e. sum_0 = c_00 + c_01 y + c_02 y^2 ...
    //     sum_1 = c_10 x + c_11 xy + c12 x y^2 + ...
    //     sum_2 = c_20 x^2 + c_21 x^2 y + c_22 x^2 y^2 + ...
    i=0;
    while (i!=degree+1){
        Expression sum_i=0;
        for (unsigned int j=0; j<param[i].size(); j++){
            sum_i = sum_i + Constant(0,1) * param[i][j] * pow(x, i) * pow(y, j);
            if (debug){
            std::cout<<"x = "<<x<<"\n";
            std::cout<<"y = "<<y<<"\n";
            std::cout<<"param[i][j] = "<<param[i][j]<<"\n";
            std::cout<<"sum_i = "<<sum_i<<"\n";
            }
        }
        i++;
        sum = sum + sum_i;
    }
    Expression amp = exp(sum);
    //amp = sum;
//    Expression amp = Constant(0, 1);
    //std::cout<<"amp = "<<amp<<"\n";
    return amp;
}