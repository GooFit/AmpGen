#ifndef AMPGEN_INTEGRATE_FP
#define AMPGEN_INTEGRATE_FP 1

#include <tuple>
#include "AmpGen/simd/utils.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen { 

  template<unsigned dim, typename fcn> std::tuple<double, double, unsigned> integrate_fp_avx2d( const fcn& F, const std::array<double,dim>& ctr, const std::array<double, dim>& wth )
  {
    static const double xl2 = 0.358568582800318073;//lambda_2
    static const double xl4 = 0.948683298050513796;//lambda_4
    static const double xl5 = 0.688247201611685289;//lambda_5
    static const double w2  = 980./6561; //weights/2^n
    static const double w4  = 200./19683;
    static const double wp2 = 245./486;//error weights/2^n
    static const double wp4 = 25./729;

    static const double wn1[14] = {     -0.193872885230909911, -0.555606360818980835,
      -0.876695625666819078, -1.15714067977442459,  -1.39694152314179743,
      -1.59609815576893754,  -1.75461057765584494,  -1.87247878880251983,
      -1.94970278920896201,  -1.98628257887517146,  -1.98221815780114818,
      -1.93750952598689219,  -1.85215668343240347,  -1.72615963013768225};

    static const double wn3[14] = {     0.0518213686937966768,  0.0314992633236803330,
      0.0111771579535639891,-0.00914494741655235473,-0.0294670527866686986,
      -0.0497891581567850424,-0.0701112635269013768, -0.0904333688970177241,
      -0.110755474267134071, -0.131077579637250419,  -0.151399685007366752,
      -0.171721790377483099, -0.192043895747599447,  -0.212366001117715794};

    static const double wn5[14] = {         0.871183254585174982e-01,  0.435591627292587508e-01,
      0.217795813646293754e-01,  0.108897906823146873e-01,  0.544489534115734364e-02,
      0.272244767057867193e-02,  0.136122383528933596e-02,  0.680611917644667955e-03,
      0.340305958822333977e-03,  0.170152979411166995e-03,  0.850764897055834977e-04,
      0.425382448527917472e-04,  0.212691224263958736e-04,  0.106345612131979372e-04};

    static const double wpn1[14] = {   -1.33196159122085045, -2.29218106995884763,
      -3.11522633744855959, -3.80109739368998611, -4.34979423868312742,
      -4.76131687242798352, -5.03566529492455417, -5.17283950617283939,
      -5.17283950617283939, -5.03566529492455417, -4.76131687242798352,
      -4.34979423868312742, -3.80109739368998611, -3.11522633744855959};

    static const double wpn3[14] = {     0.0445816186556927292, -0.0240054869684499309,
      -0.0925925925925925875, -0.161179698216735251,  -0.229766803840877915,
      -0.298353909465020564,  -0.366941015089163228,  -0.435528120713305891,
      -0.504115226337448555,  -0.572702331961591218,  -0.641289437585733882,
      -0.709876543209876532,  -0.778463648834019195,  -0.847050754458161859};

    auto eval = [&F](auto& z, const auto& p, const unsigned& ind )
    {
      auto zt = z[ind]; 
      z[ind] += p;
      auto v = F(z);
      z[ind] = zt;
      return v; 
    };
    auto eval_2 = [&F](auto& z, const auto& p, const unsigned& ind, const auto& p2, const unsigned int& ind2 )
    {
      auto zt  = z[ind]; 
      auto zt2 = z[ind2];
      z[ind]   += p;
      z[ind2]  += p2;
      auto v   = F(z);
      z[ind]   = zt;
      z[ind2]  = zt2;
      return v; 
    };


    unsigned idvaxn =0;
    std::array<float_v, dim> z; 
    double rgnvol = get_power<2,dim>::value; 
    for (unsigned j=0; j<dim; j++){
      z[j] = ctr[j]; 
      rgnvol *= wth[j]; //region volume
    }
    double sum1 = F(z).at(0);
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;

    double difmax = 0;
    static const float_v xl(   -xl4, -xl2, +xl2, +xl4);
    static const float_v dxp ( -xl4, +xl4, -xl4, +xl4);
    static const float_v dxp2( -xl4, -xl4, +xl4, +xl4);

    //loop over coordinates
    for (unsigned j=0; j != dim; j++)
    {
      auto fb  = eval(z, xl * wth[j], j ).to_array();
      auto f2 = fb[1] + fb[2];
      auto f3 = fb[0] + fb[3]; 
      sum2   += f2;//sum func eval with different weights separately
      sum3   += f3; 
      double dif     = std::abs(7*f2-f3-12*sum1);
      if (dif >= difmax) {
        difmax=dif;
        idvaxn=j;
      }
      for( unsigned k = j+1; k <dim; ++k )
        sum4 += utils::sum_elements( eval_2(z,  dxp*wth[j], j,  dxp2*wth[k], k) );
    }
    for( int j = 0 ; j < get_power<2, dim>::value; j+=4 )
    {
      for( int k = 0; k != dim; ++k ) z[k] = float_v( 
          ctr[k] + ( 0x1 & ( (j+0) >> k ) ? +1 : -1 ) * xl5*wth[k], 
          ctr[k] + ( 0x1 & ( (j+1) >> k ) ? +1 : -1 ) * xl5*wth[k], 
          ctr[k] + ( 0x1 & ( (j+2) >> k ) ? +1 : -1 ) * xl5*wth[k], 
          ctr[k] + ( 0x1 & ( (j+3) >> k ) ? +1 : -1 ) * xl5*wth[k] );
      sum5 += utils::sum_elements( F(z) );
    }

    auto rgncmp  = rgnvol*(wpn1[dim-2]*sum1+wp2*sum2+wpn3[dim-2]*sum3+wp4*sum4);
    auto rgnval  = rgnvol*(wn1[dim-2]*sum1+w2*sum2+wn3[dim-2]*sum3+w4*sum4+wn5[dim-2]*sum5);
    auto rgnerr  = std::abs(rgnval-rgncmp);//compares estim error with expected error
    return {rgnval, rgnerr, idvaxn};
  }

  template<unsigned dim, typename fcn> std::tuple<double, double, unsigned> integrate_fp_scalar( const fcn& F, const std::array<double,dim>& ctr, const std::array<double, dim>& wth )
  {
    static const double xl2 = 0.358568582800318073;//lambda_2
    static const double xl4 = 0.948683298050513796;//lambda_4
    static const double xl5 = 0.688247201611685289;//lambda_5
    static const double w2  = 980./6561; //weights/2^n
    static const double w4  = 200./19683;
    static const double wp2 = 245./486;//error weights/2^n
    static const double wp4 = 25./729;

    static const double wn1[14] = {     -0.193872885230909911, -0.555606360818980835,
      -0.876695625666819078, -1.15714067977442459,  -1.39694152314179743,
      -1.59609815576893754,  -1.75461057765584494,  -1.87247878880251983,
      -1.94970278920896201,  -1.98628257887517146,  -1.98221815780114818,
      -1.93750952598689219,  -1.85215668343240347,  -1.72615963013768225};

    static const double wn3[14] = {     0.0518213686937966768,  0.0314992633236803330,
      0.0111771579535639891,-0.00914494741655235473,-0.0294670527866686986,
      -0.0497891581567850424,-0.0701112635269013768, -0.0904333688970177241,
      -0.110755474267134071, -0.131077579637250419,  -0.151399685007366752,
      -0.171721790377483099, -0.192043895747599447,  -0.212366001117715794};

    static const double wn5[14] = {         0.871183254585174982e-01,  0.435591627292587508e-01,
      0.217795813646293754e-01,  0.108897906823146873e-01,  0.544489534115734364e-02,
      0.272244767057867193e-02,  0.136122383528933596e-02,  0.680611917644667955e-03,
      0.340305958822333977e-03,  0.170152979411166995e-03,  0.850764897055834977e-04,
      0.425382448527917472e-04,  0.212691224263958736e-04,  0.106345612131979372e-04};

    static const double wpn1[14] = {   -1.33196159122085045, -2.29218106995884763,
      -3.11522633744855959, -3.80109739368998611, -4.34979423868312742,
      -4.76131687242798352, -5.03566529492455417, -5.17283950617283939,
      -5.17283950617283939, -5.03566529492455417, -4.76131687242798352,
      -4.34979423868312742, -3.80109739368998611, -3.11522633744855959};

    static const double wpn3[14] = {     0.0445816186556927292, -0.0240054869684499309,
      -0.0925925925925925875, -0.161179698216735251,  -0.229766803840877915,
      -0.298353909465020564,  -0.366941015089163228,  -0.435528120713305891,
      -0.504115226337448555,  -0.572702331961591218,  -0.641289437585733882,
      -0.709876543209876532,  -0.778463648834019195,  -0.847050754458161859};

    auto eval = [&F](auto& z, const auto& p, const unsigned& ind )
    {
      auto zt = z[ind]; 
      z[ind] += p;
      auto v = F(z);
      z[ind] = zt;
      return v; 
    };
    auto eval_2 = [&F](auto& z, const auto& p, const unsigned& ind, const auto& p2, const unsigned int& ind2 )
    {
      auto zt  = z[ind]; 
      auto zt2 = z[ind2];
      z[ind]   += p;
      z[ind2]  += p2;
      auto v   = F(z);
      z[ind]   = zt;
      z[ind2]  = zt2;
      return v; 
    };


    unsigned idvaxn =0;
    std::array<double, dim> z = ctr;   
    double rgnvol = get_power<2,dim>::value; 
    for (unsigned j=0; j<dim; j++)  rgnvol *= wth[j]; //region volume
    double sum1 = F(z);
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;

    double difmax = 0;

    //loop over coordinates
    for (unsigned j=0; j != dim; j++)
    {
      auto f2  = eval(z, -xl2*wth[j], j) + eval(z, +xl2*wth[j], j);
      auto f3  = eval(z, -xl4*wth[j], j) + eval(z, +xl4*wth[j], j );

      sum2   += f2;//sum func eval with different weights separately
      sum3   += f3; 
      double dif     = std::abs(7*f2-f3-12*sum1);
      //storing dimension with biggest error/difference (?)
      if (dif >= difmax) {
        difmax=dif;
        idvaxn=j;
      }
      for( unsigned k = j+1; k <dim; ++k )
      {
        sum4 += eval_2(z,  xl4*wth[j],  j,  xl4*wth[k], k);
        sum4 += eval_2(z,  xl4*wth[j],  j, -xl4*wth[k], k);
        sum4 += eval_2(z,  -xl4*wth[j], j,  xl4*wth[k], k);
        sum4 += eval_2(z,  -xl4*wth[j], j, -xl4*wth[k], k);
      }
    }
    for( int j = 0 ; j != get_power<2, dim>::value; ++j )
    {
      for( int k = 0; k != dim; ++k ) z[k] = ctr[k] + ( 0x1 & ( j >> k ) ? +1 : -1 ) * xl5*wth[k];
      sum5 += F(z);
    }

    auto rgncmp  = rgnvol*(wpn1[dim-2]*sum1+wp2*sum2+wpn3[dim-2]*sum3+wp4*sum4);
    auto rgnval  = rgnvol*(wn1[dim-2]*sum1+w2*sum2+wn3[dim-2]*sum3+w4*sum4+wn5[dim-2]*sum5);
    auto rgnerr  = std::abs(rgnval-rgncmp);//compares estim error with expected error
    return {rgnval, rgnerr, idvaxn};
  }

}

#endif
