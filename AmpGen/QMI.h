#ifndef AMPGEN_QMI_H
#define AMPGEN_QMI_H

#include <memory.h>
#include <stddef.h>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Types.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Store.h"
#include "AmpGen/KeyedFunctors.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/BinDT.h"


#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <map>
#include "TMatrixD.h"
#include "TDecompSVD.h"





#ifdef _OPENMP
  #include <omp.h>
#endif
#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

namespace AmpGen{
    namespace QMI{

        std::vector<Expression> dalitz(EventType type){
        Particle mother(type.decayDescriptor(), type.finalStates());
            Tensor p0(mother.P());
            Tensor p1(mother.daughter(0)->P());
            Tensor p2(mother.daughter(1)->P());
            Tensor p3(mother.daughter(2)->P());
            Expression s12 = 
                (p1[3] + p2[3]) * (p1[3] + p2[3]) - 
                (p1[0] + p2[0]) * (p1[0] + p2[0]) - 
                (p1[1] + p2[1]) * (p1[1] + p2[1]) - 
                (p1[2] + p2[2]) * (p1[2] + p2[2]);
            Expression s13 = 
                (p1[3] + p3[3]) * (p1[3] + p3[3]) - 
                (p1[0] + p3[0]) * (p1[0] + p3[0]) - 
                (p1[1] + p3[1]) * (p1[1] + p3[1]) - 
                (p1[2] + p3[2]) * (p1[2] + p3[2]);
            Expression s23 = 
                (p2[3] + p3[3]) * (p2[3] + p3[3]) - 
                (p2[0] + p3[0]) * (p2[0] + p3[0]) - 
                (p2[1] + p3[1]) * (p2[1] + p3[1]) - 
                (p2[2] + p3[2]) * (p2[2] + p3[2]);
            
            Expression P1 = fcn::pow(p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2], 0.5);
//            Expression M1 = fcn::pow(fcn::pow(P1, 2) + fcn::pow(p1[3], 2), 0.5);
            real_t M1 = mother.daughter(0)->mass();
            real_t M2 = mother.daughter(1)->mass();
            real_t M3 = mother.daughter(2)->mass();
            Expression P2 = fcn::pow(p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2], 0.5);
 //           Expression M2 = fcn::pow(fcn::pow(P2, 2) + fcn::pow(p2[3], 2), 0.5);
            Expression P3 = fcn::pow(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2], 0.5);
//            Expression M3 = fcn::pow(fcn::pow(P3, 2) + fcn::pow(p3[3], 2), 0.5);
//            Expression M = fcn::pow(s12 + s13 + s23 - fcn::pow(M1, 2) - fcn::pow(M2, 2) - fcn::pow(M3, 2), 0.5);
            Expression cos23 = (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2])/(P2 * P3);
//            Expression theta23 = fcn::acos(cos23);
//            Expression cos23 = (fcn::pow(p1[3], 2) - fcn::pow(M1, 2) - fcn::pow(P2, 2) - fcn::pow(P3, 2))/(2 * P2 * P3);
//            Expression cos23 = (fcn::pow(M, 2) + fcn::pow(M2, 2) + fcn::pow(M3, 2) - fcn::pow(M1, 2) - 2 * M * (p2[3] + p3[3]) + 2 * p2[3] * p3[3])/(2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]));
            Expression E1 = (fcn::pow(M1, 2) + mother.mass() * mother.mass() - s23)/(2 * mother.mass());
            Expression E2 = (fcn::pow(M2, 2) + mother.mass() * mother.mass() - s13)/(2 * mother.mass());
            Expression E3 = (fcn::pow(M3, 2) + mother.mass() * mother.mass() - s12)/(2 * mother.mass());
            Expression theta23 = fcn::acos(cos23);
            Expression cos12 = (2 * E1 * E2 - (s12 - fcn::pow(M1, 2) - fcn::pow(M2, 2)))/(2 * fcn::sqrt( fcn::pow(E1, 2) - fcn::pow(M1, 2)) * fcn::sqrt(fcn::pow(E2, 2) - fcn::pow(M2, 2)));
            Expression cos13 = (2 * E1 * E3 - (s13 - fcn::pow(M1, 2) - fcn::pow(M3, 2)))/(2 * fcn::sqrt( fcn::pow(E1, 2) - fcn::pow(M1, 2)) * fcn::sqrt(fcn::pow(E3, 2) - fcn::pow(M3, 2)));
            //Expression cos13 = (2 * E1 * E3 - (s12 - fcn::pow(M1, 2) - fcn::pow(M3, 2)))/(2 * fcn::sqrt(P1 * P3));
            Expression theta = (fcn::acos(cos12) - fcn::acos(cos13));
            return {s12, s13, s23, theta};

        }



        std::vector<Expression> rotateDalitz(EventType& type){
            //std::vector<Expression> Phi(dalitz(type));
            real_t c1 = -3.1171885586526695;
            real_t m1 = 2.23407421671132946;
            real_t c2 = -9.54231895051727e-05;
            real_t m2 = 0.8051636393861085;
             
            std::vector<Expression> Phi(dalitz(type));
            Expression z1 = (Phi[1] + Phi[0])/2;
            Expression z2 = (Phi[1] - Phi[0])/2;
            Expression w1 = m1 * z1 + c1;
            Expression w2 = m2 * z2 + c2;
            if (NamedParameter<bool>("PhaseCorrection::stretchAntiSym", false)){
                    w2 = Expression(NamedParameter<real_t>("PhaseCorrection::stretchAntiSym_A", 2)) * w2/(w1 + 1 - Expression(NamedParameter<real_t>("PhaseCorrection::stretchAntiSym_epsilon", 0.01)));
            }
            if (NamedParameter<bool>("PhaseCorrection::useSquareDalitz", false)){
                    z1 = Phi[2];
                    z2 = Phi[3];
                    real_t z1Min = std::pow(2 * 0.13957039, 2);
                    real_t z1Max = std::pow(1.86484 - 0.497611, 2);
//                    real_t z2Min = -M_PI;
//                    real_t z2Max = M_PI;
//                    m1 = 2/(z1Max - z1Min);
//                    c1 = (z1Max + z1Min)/(z1Min - z1Max);
//                    m2 = 2/(z2Max - z2Min);
//                    c2 = (z2Max + z2Min)/(z2Min - z2Max);
//                    w1 = m1 * z1 + c1;
//                    w1 = (1/M_PI) * fcn::acos(2 * (z1 - z1Min)/(z1Max - z1Min) - 1);                    
                    w1 = 2 * (z1 - z1Min)/(z1Max - z1Min) - 1;
                    w2 = z2/M_PI;//fcn::sin(z2);//m2 * z2 + c2;
                    }
            return {w1, w2};



        }

        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        std::vector<ce> xy(EventType type, MinuitParameterSet& MPS){
            std::vector<Expression> w(rotateDalitz(type));
            //CompiledExpression<real_t, real_t*, real_t*>comp_w1(w1, "w1", type.getEventFormat(), {} ,&MPS);
            ce comp_w1(w[0], "w1", &MPS, type.getEventFormat());
            //CompiledExpression<real_t, real_t*, real_t*>comp_w2(w2, "w2", type.getEventFormat(), {},&MPS);
            ce comp_w2(w[1], "w2", &MPS, type.getEventFormat());
            comp_w1.prepare(); comp_w1.compile();
            comp_w2.prepare(); comp_w2.compile();
            return {comp_w1, comp_w2};
                    
        }

        Expression legendre(Expression& x, size_t n){
            if (n==0) return 1;
            if (n==1) return x;
            Expression num = (2 - 1/n) * x * legendre(x, n-1) - (1 - 2/n) * legendre(x, n-2);

            return num;
        }

        //template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        Expression Pij(EventType& type, size_t i, size_t j, std::string polyType){
            std::vector<Expression> w(rotateDalitz(type));
            Expression Pi, Pj;
            
            if(polyType == "antiSym_legendre"){
                Pi = legendre(w[0], i);
                Pj = legendre(w[1], j);
            }
            
            if (polyType=="antiSym_simple"){
                Pi = fcn::pow(w[0], i);
                Pj = fcn::pow(w[1], j);
            }
 
            return Pi * Pj;
        }
//        Expression polyGaussian(EventType type){
//            std::vector<Expression> Phi(dalitz(type));
//            std::vector<Expression> w(rotateDalitz(type));
//            
 //       }


        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        ce cPij(EventType& type, size_t i, size_t j, MinuitParameterSet& MPS, std::string polyType){
            Expression pij(Pij(type, i, j, polyType));
            ce cpij(pij, "PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(j), type.getEventFormat(), &MPS);
            cpij.prepare(); cpij.compile();
            return cpij;
        }


        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        std::map<std::string, ce> cPhaseCorrection(EventType& type, MinuitParameterSet& MPS){
            size_t O(NamedParameter<size_t>("PhaseCorrection::Order", 4));
            std::string polyType(NamedParameter<std::string>("PhaseCorrection::PolyType", "antiSym_legendre"));
            std::map<std::string, ce> m;
            if (polyType=="Gaussian"){
                INFO("Doing a Gaussian Bias so no compiling needed");


            }
            if (polyType=="Gaussian2"){
                INFO("Doing 2x Gaussian Bias so no compiling needed");


            }
 
            else{ 
                for (size_t i=0;i<O+1;i++){
//                   real_t j_max = 0.5 * (real_t)(O + 1 - i) - 0.5;
//                    size_t j_max_int = std::floor(j_max) + 1;
//                    INFO("jmax = "<<j_max_int);
                    for (size_t j=0;j<O +1 - i;j++){
    //                    std::string ij = std::to_string(i) + "_" + std::to_string(2 * j + 1);
                        if (i + 2 * j + 1 <= O){                   
                            INFO("Compiling O"<<i<<"_"<<2 * j + 1);
                            ce cpij(cPij(type, i, 2 * j + 1, MPS, polyType));
                            std::string ij = cpij.name();
                            std::pair<std::string, ce> p({ij,cpij});
                            m.insert(p);
                        }
                    }
                }
            }
            return m;
        }

   /*
        CompiledExpression<void, real_t, real_t> x(EventType sigType){
            return ce;
        }
        CompiledExpression<void, real_t, real_t> y(EventType sigType){
            return ce;
        }
    */ 

        real_t Cij(size_t i, size_t j, MinuitParameterSet& MPS){
            return MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(j)]->mean();
        }

        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >

        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t phaseCorrection(const Event& evt, ce& x, ce& y, MinuitParameterSet& MPS){
            size_t O(NamedParameter<size_t>("PhaseCorrection::Order", 0));
        //    std::string polyType(NamedParameter<std::string>("PhaseCorrection", "antiSym_legendre"));
            real_t f = 0;
    //        real_t * a(evt);
    //        real_t X(x({}, evt.address()));

            for (size_t i=0;i<O + 1;i++){
                for (size_t j=0;j<O + 1 - i; j++){
                    //f += MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j + 1)]->mean() * std::legendre(i, x({}, evt.address())) * std::legendre(2 * j + 1, y({}, evt.address())); 
                    //f += MPS["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j + 1)]->mean() ;//* std::legendre(i, x({}, evt.address())) * std::legendre(2 * j + 1, y({}, evt.address())); 
                    f += Cij(i, 2 * j + 1, MPS) * std::legendre(i, x(evt.address())) * std::legendre(2 * j + 1, y(evt.address()));
                }
            }
            return f;
        }


        real_t GaussDalitz(const Event& evt, MinuitParameterSet& MPS){
//            real_t rt(0);
//            real_t A(MPS["PhaseCorrection::Gauss::A"]->mean());
//            real_t w_0(MPS["PhaseCorrection::Gauss::wminus_0"]->mean());
//            real_t muplus(MPS["PhaseCorrection::Gauss::mu_+"]->mean());
//            real_t sigmaplus(MPS["PhaseCorrection::Gauss::sigma_+"]->mean());
//            real_t muminus(MPS["PhaseCorrection::Gauss::mu_-"]->mean());
//            real_t sigmaminus(MPS["PhaseCorrection::Gauss::sigma_-"]->mean());
//            real_t s01(evt.s(0, 1));
//            real_t s02(evt.s(0, 2));
//            real_t c2 = -9.54231895051727e-05;
//            real_t m2 = 0.8051636393861085;
//            real_t z2 = (s01 - s02)/2;
//            real_t w2 = m2 * z2 + c2;
//            if (s01>s02) {
//                rt =  A * std::erf(w2/w_0) * std::exp(-std::pow(s01 - muplus, 2)/(2 * sigmaplus)) * std::exp(-std::pow(s02 - muminus, 2)/(2 * sigmaminus));
//            }
//            else{
//                rt =  A * std::erf(w2/w_0) * std::exp(-std::pow(s02 - muplus, 2)/(2 * sigmaplus)) * std::exp(-std::pow(s01 - muminus, 2)/(2 * sigmaminus));
//
//            }
//           // INFO("Gauss = "<<rt);
//            return rt; 

      real_t A(MPS["PhaseCorrection::Gauss::A"]->mean());
        real_t w_0(MPS["PhaseCorrection::Gauss::wminus_0"]->mean());
        real_t muplus(MPS["PhaseCorrection::Gauss::mu_+"]->mean());
        real_t sigmaplus(MPS["PhaseCorrection::Gauss::sigma_+"]->mean());
        real_t muminus(MPS["PhaseCorrection::Gauss::mu_-"]->mean());
        real_t sigmaminus(MPS["PhaseCorrection::Gauss::sigma_-"]->mean());
        real_t s01(evt.s(0, 1));
        real_t s02(evt.s(0, 2));
        real_t c2 = -9.54231895051727e-05;
        real_t m2 = 0.8051636393861085;
        real_t z2 = (s01 - s02)/2;
        real_t w2 = m2 * z2 + c2;


        real_t gauss_pp = std::pow(s01 - muplus, 2);
        real_t gauss_pm = std::pow(s01 - muminus, 2);
        real_t gauss_mm = std::pow(s02 - muminus, 2);
        real_t gauss_mp = std::pow(s02 - muplus, 2);
        real_t z_pp = gauss_pp/(2 * std::pow(sigmaplus, 2));
        real_t z_pm = gauss_pm/(2 * std::pow(sigmaminus, 2));
        real_t z_mp = gauss_mp/(2 * std::pow(sigmaplus, 2));
        real_t z_mm = gauss_mm/(2 * std::pow(sigmaminus,2));
        real_t exp_pp = std::exp(-z_pp);
        real_t exp_pm = std::exp(-z_pm);
        real_t exp_mp = std::exp(-z_mp);
        real_t exp_mm = std::exp(-z_mm);


        real_t g = 0;
        if (s01>s02) g = A*std::erf(w2/w_0)*exp_pp * exp_mm;
//        if (s01>s02) g = A*exp_pp * exp_mm;
        if (s01<s02) g = A*std::erf(w2/w_0)*exp_pm * exp_mp;
       // if (s01<s02) g = -A*exp_pm * exp_mp;
        return g;
        }
        
        real_t Gauss2Dalitz(const Event& evt, MinuitParameterSet& MPS){
//            real_t rt(0);
//            real_t A(MPS["PhaseCorrection::Gauss::A"]->mean());
//            real_t w_0(MPS["PhaseCorrection::Gauss::wminus_0"]->mean());
//            real_t muplus(MPS["PhaseCorrection::Gauss::mu_+"]->mean());
//            real_t sigmaplus(MPS["PhaseCorrection::Gauss::sigma_+"]->mean());
//            real_t muminus(MPS["PhaseCorrection::Gauss::mu_-"]->mean());
//            real_t sigmaminus(MPS["PhaseCorrection::Gauss::sigma_-"]->mean());
//            real_t s01(evt.s(0, 1));
//            real_t s02(evt.s(0, 2));
//            real_t c2 = -9.54231895051727e-05;
//            real_t m2 = 0.8051636393861085;
//            real_t z2 = (s01 - s02)/2;
//            real_t w2 = m2 * z2 + c2;
//            if (s01>s02) {
//                rt =  A * std::erf(w2/w_0) * std::exp(-std::pow(s01 - muplus, 2)/(2 * sigmaplus)) * std::exp(-std::pow(s02 - muminus, 2)/(2 * sigmaminus));
//            }
//            else{
//                rt =  A * std::erf(w2/w_0) * std::exp(-std::pow(s02 - muplus, 2)/(2 * sigmaplus)) * std::exp(-std::pow(s01 - muminus, 2)/(2 * sigmaminus));
//
//            }
//           // INFO("Gauss = "<<rt);
//            return rt; 

        real_t A1(MPS["PhaseCorrection::Gauss2::A1"]->mean());
        real_t A2(MPS["PhaseCorrection::Gauss2::A2"]->mean());
        real_t w_01(MPS["PhaseCorrection::Gauss2::wminus1_0"]->mean());
        real_t w_02(MPS["PhaseCorrection::Gauss2::wminus2_0"]->mean());
        real_t muplus1(MPS["PhaseCorrection::Gauss2::mu1_+"]->mean());
        real_t muplus2(MPS["PhaseCorrection::Gauss2::mu2_+"]->mean());
        real_t sigmaplus1(MPS["PhaseCorrection::Gauss2::sigma1_+"]->mean());
        real_t sigmaplus2(MPS["PhaseCorrection::Gauss2::sigma2_+"]->mean());
        real_t muminus1(MPS["PhaseCorrection::Gauss2::mu1_-"]->mean());
        real_t muminus2(MPS["PhaseCorrection::Gauss2::mu2_-"]->mean());
        real_t sigmaminus1(MPS["PhaseCorrection::Gauss2::sigma1_-"]->mean());
        real_t sigmaminus2(MPS["PhaseCorrection::Gauss2::sigma2_-"]->mean());
        real_t s01(evt.s(0, 1));
        real_t s02(evt.s(0, 2));
        real_t c2 = -9.54231895051727e-05;
        real_t m2 = 0.8051636393861085;
        real_t z2 = (s01 - s02)/2;
        real_t w2 = m2 * z2 + c2;


        real_t gauss_pp1 = std::pow(s01 - muplus1, 2);
        real_t gauss_pm1 = std::pow(s01 - muminus1, 2);
        real_t gauss_mm1 = std::pow(s02 - muminus1, 2);
        real_t gauss_mp1 = std::pow(s02 - muplus1, 2);
        real_t z_pp1 = gauss_pp1/(2 * std::pow(sigmaplus1, 2));
        real_t z_pm1 = gauss_pm1/(2 * std::pow(sigmaminus1, 2));
        real_t z_mp1 = gauss_mp1/(2 * std::pow(sigmaplus1, 2));
        real_t z_mm1 = gauss_mm1/(2 * std::pow(sigmaminus1,2));
        real_t exp_pp1 = std::exp(-z_pp1);
        real_t exp_pm1 = std::exp(-z_pm1);
        real_t exp_mp1 = std::exp(-z_mp1);
        real_t exp_mm1 = std::exp(-z_mm1);

        real_t gauss_pp2 = std::pow(s01 - muplus2, 2);
        real_t gauss_pm2 = std::pow(s01 - muminus2, 2);
        real_t gauss_mm2 = std::pow(s02 - muminus2, 2);
        real_t gauss_mp2 = std::pow(s02 - muplus2, 2);
        real_t z_pp2 = gauss_pp2/(2 * std::pow(sigmaplus2, 2));
        real_t z_pm2 = gauss_pm2/(2 * std::pow(sigmaminus2, 2));
        real_t z_mp2 = gauss_mp2/(2 * std::pow(sigmaplus2, 2));
        real_t z_mm2 = gauss_mm2/(2 * std::pow(sigmaminus2,2));
        real_t exp_pp2 = std::exp(-z_pp2);
        real_t exp_pm2 = std::exp(-z_pm2);
        real_t exp_mp2 = std::exp(-z_mp2);
        real_t exp_mm2 = std::exp(-z_mm2);





        real_t g1 = 0;
        real_t g2 = 0;
        if (s01>s02) g1 = A1*std::erf(w2/w_01)*exp_pp1* exp_mm1;
        if (s01>s02) g2 = A2*std::erf(w2/w_02)*exp_pp2* exp_mm2;
//        if (s01>s02) g = A*exp_pp * exp_mm;
        if (s01<s02) g1 = A1*std::erf(w2/w_01)*exp_pm1 * exp_mp1;
        if (s01<s02) g2 = A2*std::erf(w2/w_02)*exp_pm2 * exp_mp2;
       // if (s01<s02) g = -A*exp_pm * exp_mp;
        return g1 + g2;
        }
 

        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t phaseCorrection(const Event& evt, std::map<std::string, ce>& cPij, MinuitParameterSet& MPS){
            real_t rt=0;
            if (NamedParameter<std::string>("PhaseCorrection::PolyType", "antiSym_legendre")  == "Gaussian" ){
               rt = GaussDalitz(evt, MPS);
            }
            else if (NamedParameter<std::string>("PhaseCorrection::PolyType", "antiSym_legendre")  == "Gaussian2" ){
               rt = Gauss2Dalitz(evt, MPS);
            }
            else{
                for (auto& p : cPij) rt += MPS[p.first]->mean() * p.second(evt.address());
            }
           return rt; 
            //return std::asin(std::sin(rt));
       }




        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        complex_t AC(CoherentSum& A, CoherentSum& Abar, EventList& evts, ce& x, ce& y, MinuitParameterSet& MPS){
            complex_t rt =0;
            for (auto& evt : evts){
                complex_t a(A.getValNoCache(evt) );
                complex_t abar(Abar.getValNoCache(evt));
                real_t rD = std::abs(abar/a);
                real_t dd = std::arg(abar/a);
                real_t f = phaseCorrection(evt, x, y, MPS);
                rt += std::abs(a) * std::abs(abar) * complex_t(cos(dd + f), sin(dd + f));
            }
            return rt;
        }

        real_t dd(const Event& evt,CoherentSum& A, CoherentSum& Abar){
                //complex_t a(A.getValNoCache(evt) );
                //complex_t abar(Abar.getValNoCache(evt));
                //real_t rD = std::abs(abar/a);
//                real_t _dd = std::arg(abar) - std::arg(a);
            return std::arg(A.getValNoCache(evt) * std::conj(Abar.getValNoCache(evt))) ;
            //return -M_PI +  std::arg(A.getValNoCache(evt) * std::conj(Abar.getValNoCache(evt))) ;
//            return std::arg(std::conj(A.getValNoCache(evt)) * Abar.getValNoCache(evt)) ;
        }
        std::map<std::string, std::vector<real_t> > AmpArrays(const EventList_type& mcSig, CoherentSum& A, CoherentSum& Abar){
           
            std::vector<real_t> dds(mcSig.size());
            std::vector<real_t> absA(mcSig.size());
            std::vector<real_t> absAbar(mcSig.size());
            #pragma omp parallel for
            for (unsigned i=0;i<mcSig.size();++i){
                absA[i] = std::abs(A.getValNoCache(mcSig[i]));
                absAbar[i] = std::abs(Abar.getValNoCache(mcSig[i]));
                dds[i] = dd(mcSig[i], A, Abar);
            }
            
            return std::map<std::string, std::vector<real_t> >({{"dd", dds}, {"A", absA}, {"Abar", absAbar}} );
        }


 

        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t cosTerm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, ce& x, ce& y, EventList_type& mcSig){ 
            real_t C=0;
            //#pragma omp parallel for reduction(+:C)
            for (auto& evt : mcSig) C +=  std::abs(A.getValNoCache(evt) * Abar.getValNoCache(evt)) * cos(dd(evt, A, Abar) + phaseCorrection(evt, x, y, MPS));
            
            return C;
        }
        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t cosTerm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig){ 
            real_t C=0;
            
           #pragma omp parallel for reduction(+:C) 
           for(unsigned i=0;i<mcSig.size();++i){
                C += std::abs(A.getValNoCache(mcSig[i]) * Abar.getValNoCache(mcSig[i])) * cos(dd(mcSig[i], A, Abar) + phaseCorrection(mcSig[i], cPij, MPS));
            }
            
            return C/(real_t)mcSig.size();
        }

        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t cosTerm(std::map<std::string, std::vector<real_t> >& amps, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig){ 
            real_t C=0;
         //   INFO("Starting CosTerm Function, before OMP loop");
         //   INFO("Have "<<mcSig.size()<<" MC events for Cos Term"); 
         //   INFO("Have "<<amps["A"].size()<<" saved |A|'s from MC");
            
           #pragma omp parallel for reduction(+:C) 
           for(unsigned i=0;i<mcSig.size();++i){
                C += amps["A"][i] * amps["Abar"][i] * cos(amps["dd"][i] + phaseCorrection(mcSig[i], cPij, MPS));
           } 
            return C/(real_t)mcSig.size();
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t cosTerm(std::map<std::string, std::vector<real_t> >* amps, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig){ 
            real_t C=0;
           // INFO("Starting CosTerm Function, before OMP loop");
           // INFO("Have "<<mcSig.size()<<" MC events for Cos Term"); 
           // INFO("Have "<<(*amps)["A"].size()<<" saved |A|'s from MC");
            
           #pragma omp parallel for reduction(+:C) 
           for(unsigned i=0;i<mcSig.size();++i){
                C += (*amps)["A"][i] * (*amps)["Abar"][i] * cos((*amps)["dd"][i] + phaseCorrection(mcSig[i], cPij, MPS));
           } 
            return C/(real_t)mcSig.size();
        }
 
 
        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t sinTerm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, ce& x, ce& y, EventList_type& mcSig){ 
            real_t S=0;
            
            //#pragma omp parallel for reduction(+:S) 
            for (auto& evt : mcSig) S +=  std::abs(A.getValNoCache(evt) * Abar.getValNoCache(evt)) * sin(dd(evt, A, Abar) + phaseCorrection(evt, x, y, MPS));
            
            return S;
        }
        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t sinTerm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig){ 
            real_t S=0;
            
            #pragma omp parallel for reduction(+:S) 
            for (unsigned i=0;i<mcSig.size();++i) S += std::abs(A.getValNoCache(mcSig[i]) * Abar.getValNoCache(mcSig[i])) * sin(dd(mcSig[i], A, Abar) + phaseCorrection(mcSig[i], cPij, MPS));
            
            return S/(real_t)mcSig.size();
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t sinTerm(std::map<std::string, std::vector<real_t> >& amps, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig){ 
            real_t S=0;
            
           #pragma omp parallel for reduction(+:S) 
           for(unsigned i=0;i<mcSig.size();++i){
                S += amps["A"][i] * amps["Abar"][i] * sin(amps["dd"][i] + phaseCorrection(mcSig[i], cPij, MPS));
           } 
            return S/(real_t)mcSig.size();
        }

        real_t cosTermNoCorrection(std::map<std::string, std::vector<real_t> >& amps, EventList_type& mcSig){
            real_t C =0;
           #pragma omp parallel for reduction(+:C) 
           for(unsigned i=0;i<mcSig.size();++i){
                C += amps["A"][i] * amps["Abar"][i] * cos(amps["dd"][i]);
           } 
            return C/(real_t)mcSig.size();
 
        }

        real_t sinTermNoCorrection(std::map<std::string, std::vector<real_t> >& amps, EventList_type& mcSig){
            real_t S =0;
           #pragma omp parallel for reduction(+:S) 
           for(unsigned i=0;i<mcSig.size();++i){
                S += amps["A"][i] * amps["Abar"][i] * sin(amps["dd"][i]);
           } 
            return S/(real_t)mcSig.size();
 
        }







        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t corr_norm(CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, EventList_type& mcSig, ce& x, ce& y, bool sameTag){
            real_t pureTerms = A.norm() * Bbar.norm() + Abar.norm() * B.norm();
            real_t C = cosTerm(A, Abar, MPS, x, y, mcSig);
            if (sameTag){
                real_t S = sinTerm(A, Abar, MPS, x, y, mcSig);
                return 2 * A.norm() * Abar.norm() - 2 * std::pow(C, 2) - 2 * std::pow(C, 2);
            }
            else{
                return A.norm() * Bbar.norm() + Abar.norm() * B.norm() - 2 * C;
            }   
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t corr_norm(CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, EventList_type& mcSig, std::map<std::string, ce>& cPij, bool sameTag){
            real_t pureTerms = A.norm() * Bbar.norm() + Abar.norm() * B.norm();
            real_t C = cosTerm(A, Abar, MPS, cPij, mcSig);
            if (sameTag){
                real_t S = sinTerm(A, Abar, MPS, cPij, mcSig);
                return 2 * A.norm() * Abar.norm() - 2 * std::pow(C, 2) - 2 * std::pow(S, 2);
            }
            else{
                return pureTerms - 2 * C;
            }   
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t corr_norm(CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, EventList_type& mcSig, std::map<std::string, ce>& cPij, bool sameTag, std::map<std::string, std::vector<real_t> >& amps){
            real_t pureTerms = A.norm() * Bbar.norm() + Abar.norm() * B.norm();
            real_t C = cosTerm(amps, MPS, cPij, mcSig);
            if (sameTag){
                real_t S = sinTerm(amps, MPS, cPij, mcSig);
                return 2 * A.norm() * Abar.norm() - 2 * std::pow(C, 2) - 2 * std::pow(S, 2);
            }
            else{
                return pureTerms - 2 * C;
            }   
        }
        

        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t corr_norm(CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, EventList_type& mcSig, EventList_type& mcTag, std::map<std::string, ce>& cPij, bool sameTag, 
                        std::map<std::string, std::vector<real_t> >& ampsSig, std::map<std::string, std::vector<real_t> >& ampsTag){
    
            real_t pureTerms = A.norm() * Bbar.norm() + Abar.norm() * B.norm();
            real_t C = cosTerm(ampsSig, MPS, cPij, mcSig);
            real_t CTag = cosTermNoCorrection(ampsTag, mcTag);
            real_t S = sinTerm(ampsSig, MPS, cPij, mcSig);
            real_t STag = sinTermNoCorrection(ampsTag, mcTag);
            if (sameTag){
                return 2 * A.norm() * Abar.norm() - 2 *C * C - 2 * S *S ;
            }
            else{
                return pureTerms - 2* C * CTag - 2 * S * STag;
            }
        }





        complex_t ckm_zB(MinuitParameterSet& MPS, int gammaSign){
            bool isCart(NamedParameter<bool>("CKM::Cartesian", true));
            if (isCart){
                if (gammaSign>0){
                    return complex_t(MPS["CKM::x+"]->mean(), MPS["CKM::y+"]->mean());
                }
                else{
                    return complex_t(MPS["CKM::x-"]->mean(), MPS["CKM::y-"]->mean());
                }
            }
            else{
                real_t angleB(MPS["CKM::deltaB"]->mean() + gammaSign *  MPS["CKM::gamma"]->mean());
                return MPS["CKM::rB"]->mean() * complex_t(cos(angleB), sin(angleB));
            }
        }


        //template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t ckm_norm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, EventList_type& mcSig, ce& x, ce& y, int gammaSign){  
            real_t C = cosTerm(A, Abar, MPS, x, y, mcSig);
            real_t S = sinTerm(A, Abar, MPS, x, y, mcSig);
            complex_t zB(ckm_zB(MPS, gammaSign));
            if (gammaSign>0){ 
                return Abar.norm() + A.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
            else{
                return A.norm() + Abar.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t ckm_norm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, EventList_type& mcSig, std::map<std::string, ce>& cPij, int gammaSign){  
            real_t C = cosTerm(A, Abar, MPS, cPij, mcSig);
            real_t S = sinTerm(A, Abar, MPS, cPij, mcSig);
            complex_t zB(ckm_zB(MPS, gammaSign));
            if (gammaSign>0){ 
                return Abar.norm() + A.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
            else{
                return A.norm() + Abar.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        real_t ckm_norm(CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, EventList_type& mcSig, std::map<std::string, ce>& cPij, int gammaSign, std::map<std::string, std::vector<real_t> >& amps){  
            real_t C = cosTerm(amps, MPS, cPij, mcSig);
            real_t S = sinTerm(amps, MPS, cPij, mcSig);
            complex_t zB(ckm_zB(MPS, gammaSign));
            if (gammaSign>0){ 
                return Abar.norm() + A.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
            else{
                return A.norm() + Abar.norm() * std::norm(zB) + 2 * (std::real(zB)*C + gammaSign * std::imag(zB) * S);
            }
        }



        real_t probCorr_unnorm_nonInt(const Event& evt1, const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, bool sameTag){
            if (sameTag){
                return std::norm(A.getValNoCache(evt1) * Abar.getValNoCache(evt2)) + std::norm(Abar.getValNoCache(evt1) * A.getValNoCache(evt2));
            }
            else{
                return std::norm(A.getValNoCache(evt1) * Bbar.getValNoCache(evt2)) + std::norm(Abar.getValNoCache(evt1) * B.getValNoCache(evt2));
            }
        }
        real_t probCorr_unnorm_nonInt(size_t idx, std::map<std::string, std::vector<real_t> >& ampsA, std::map<std::string, std::vector<real_t> >& ampsB, bool sameTag){
            return std::norm(ampsA["A"][idx] * ampsB["Abar"][idx]) + std::norm(ampsA["Abar"][idx] * ampsB["A"][idx]);
        }

        real_t probCKM_unnorm_nonInt(const Event& evt, CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, int gammaSign){
            complex_t zB(ckm_zB(MPS, gammaSign));
            if (gammaSign>0){
                return std::norm(Abar.getValNoCache(evt)) + std::norm(A.getValNoCache(evt)) * std::norm(zB);
            }
            else{
                return std::norm(A.getValNoCache(evt)) + std::norm(Abar.getValNoCache(evt)) * std::norm(zB);
            }
        }
        real_t probCKM_unnorm_nonInt(size_t idx, std::map<std::string, std::vector<real_t> >& amps, MinuitParameterSet& MPS, int gammaSign){
            complex_t zB(ckm_zB(MPS, gammaSign));
            if (gammaSign>0){
                return std::norm(amps["Abar"][idx]) + std::norm(amps["A"][idx] * zB);
            }
            else{
                return std::norm(amps["A"][idx]) + std::norm(amps["Abar"][idx] * zB);
            }
        }

        



        
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCorr_unnorm_Int(const Event& evt1, const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, ce& x, ce& y, bool sameTag){
            real_t correction1 = phaseCorrection(evt1, x, y, MPS);

            if (sameTag){

                real_t correction2 = phaseCorrection(evt2, x, y, MPS);
                return - 2 * std::abs(A.getValNoCache(evt1) * Abar.getValNoCache(evt2) * A.getValNoCache(evt2) * Abar.getValNoCache(evt1)) * cos(dd(evt1, A, Abar) - dd(evt2, A, Abar) + correction1 - correction2);
            }
            else{
                return - 2 * std::abs(A.getValNoCache(evt1) * Bbar.getValNoCache(evt2) * B.getValNoCache(evt2) * Abar.getValNoCache(evt1)) * cos(dd(evt1, A, Abar) - dd(evt2, B, Bbar) + correction1);
            }
        }
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCorr_unnorm_Int(const Event& evt1,  const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, bool sameTag){
            real_t correction1 = phaseCorrection(evt1, cPij, MPS);

            if (sameTag){

                real_t correction2 = phaseCorrection(evt2, cPij, MPS);
                return - 2 * std::abs(A.getValNoCache(evt1) * Abar.getValNoCache(evt2) * A.getValNoCache(evt2) * Abar.getValNoCache(evt1)) * cos(dd(evt1, A, Abar) - dd(evt2, A, Abar) + correction1 - correction2);
            }
            else{
                return - 2 * std::abs(A.getValNoCache(evt1) * Bbar.getValNoCache(evt2) * B.getValNoCache(evt2) * Abar.getValNoCache(evt1)) * cos(dd(evt1, A, Abar) - dd(evt2, B, Bbar) + correction1);
            }
        }
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCorr_unnorm_Int(EventList& sigList, EventList& tagList, size_t idx, std::map<std::string, std::vector<real_t> >& ampsA, std::map<std::string, std::vector<real_t> >& ampsB, MinuitParameterSet& MPS, 
                                    std::map<std::string, ce>& cPij, bool sameTag){
            real_t correction1 = phaseCorrection(sigList[idx], cPij, MPS);

            if (sameTag){

                real_t correction2 = phaseCorrection(tagList[idx], cPij, MPS);
                return (-2 * ampsA["A"][idx] * ampsA["Abar"][idx] * ampsB["A"][idx] * ampsB["Abar"][idx] * cos(ampsA["dd"][idx] - ampsB["dd"][idx] + correction1 - correction2));
            }
            else{
                //return - 2 * std::abs(A.getValNoCache(evt1) * Bbar.getValNoCache(evt2) * B.getValNoCache(evt2) * Abar.getValNoCache(evt1)) * cos(dd(evt1, A, Abar) - dd(evt2, B, Bbar) + correction1);
                return (-2 * ampsA["A"][idx] * ampsA["Abar"][idx] * ampsB["A"][idx] * ampsB["Abar"][idx] * cos(ampsA["dd"][idx] - ampsB["dd"][idx] + correction1));
            }
        }



        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCKM_unnorm_Int(const Event& evt, CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, ce& x, ce& y, int gammaSign){
            complex_t zB(ckm_zB(MPS, gammaSign));
            real_t correction = phaseCorrection(evt, x, y, MPS);
            return 2 * std::abs(A.getValNoCache(evt) * Abar.getValNoCache(evt)) *(std::real(zB) *  cos(dd(evt, A, Abar) + correction) - gammaSign * std::imag(zB) * sin(dd(evt, A, Abar) + correction));

        }
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCKM_unnorm_Int(const Event& evt, CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, int gammaSign){
            complex_t zB(ckm_zB(MPS, gammaSign));
            real_t correction = phaseCorrection(evt, cPij, MPS);
            //return 2 * std::abs(A.getValNoCache(evt) * Abar.getValNoCache(evt)) *(std::real(zB) *  cos(dd(evt, A, Abar) + correction) - gammaSign * std::imag(zB) * sin(dd(evt, A, Abar) + correction));
            return 2 * std::abs(zB) * std::abs(A.getValNoCache(evt) * Abar.getValNoCache(evt)) * cos(dd(evt, A, Abar) + correction - gammaSign *  std::arg(zB));

        }

        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCKM_unnorm_Int(EventList& list, size_t idx, std::map<std::string, std::vector<real_t> >& amps, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, int gammaSign){
            complex_t zB(ckm_zB(MPS, gammaSign));
            real_t correction = phaseCorrection(list[idx], cPij, MPS);
            return 2 * std::abs(zB) * amps["A"][idx] * amps["Abar"][idx] * cos(amps["dd"][idx] + correction - gammaSign * std::arg(zB));
        }



        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCorr_unnorm(const Event& evt1, const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, ce& x, ce& y, bool sameTag){
            return probCorr_unnorm_nonInt(evt1, evt2, A, Abar, B, Bbar, sameTag) + probCorr_unnorm_Int(evt1, evt2, A, Abar, B, Bbar, MPS, x, y, sameTag);
        }
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCorr_unnorm(const Event& evt1, const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, bool sameTag){
            return probCorr_unnorm_nonInt(evt1, evt2, A, Abar, B, Bbar, sameTag) + probCorr_unnorm_Int(evt1, evt2, A, Abar, B, Bbar, MPS, cPij, sameTag);
        }
 
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCKM_unnorm(const Event& evt, CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, ce& x, ce& y, int gammaSign){
            return probCKM_unnorm_nonInt(evt, A, Abar, MPS, gammaSign) + probCKM_unnorm_Int(evt, A, Abar, MPS, x, y, gammaSign);
        }
        template <class ce=CompiledExpression<real_t,real_t*,real_t*> >
        real_t probCKM_unnorm(const Event& evt, CoherentSum& A, CoherentSum& Abar, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, int gammaSign){
            return probCKM_unnorm_nonInt(evt, A, Abar, MPS, gammaSign) + probCKM_unnorm_Int(evt, A, Abar, MPS, cPij, gammaSign);
        }
        double getMax(const EventList& events)
        {
            double max = 0.;
            #pragma omp parallel for reduction (max:max)
            //for ( const auto& evt : events ) 
            for(unsigned i=0;i<events.size();++i)
            {
                auto value             = events[i].genPdf() ; // pdf(evt) / evt.genPdf();
                //INFO("value = "<<value<<" max = "<<max);
                if ( value > max ) max = value;
            }
            INFO( "Returning normalisation constant = " << max ); 
            return max;
        }

        template <typename eventList_t, typename pdf_t> void generatePsi3770( pdf_t& pdf, eventList_t& list1, eventList_t& list2, const size_t& N, size_t seed, size_t generatorBlock )
        {    
            TRandom3 rndm(seed);
            //  double maxProb   = m_normalise ? 0 : 1;
            double maxProb = 1;
            auto size0       = list1.size();
            double totalGenerated = 0; 
            Generator<> g1(list1.eventType());
            Generator<> g2(list2.eventType());
            //pdf.reset( true );
            ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
            ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
            std::vector<bool> efficiencyReport(generatorBlock,false); 
            EventList mc1( list1.eventType() );
            EventList mc2( list2.eventType() );
            INFO("Start generation of "<<N<<" "<<mc1.eventType()<<" vs "<<mc2.eventType()<<" events");
            while ( list1.size() - size0 < N ) {
                t_phsp.start();
                g1.fillEventListPhaseSpace(mc1, generatorBlock);
                g2.fillEventListPhaseSpace(mc2, generatorBlock);
                t_phsp.stop();
                t_eval.start();
            //    pdf.setEvents( mc1 );
            //    pdf.prepare();
                auto previousSize = list1.size();

                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for ( size_t block=0; block < mc1.nBlocks(); ++block )
                { 
                //mc1.setGenPDF(block, pdf(mc1.block(block), block) / mc1.genPDF(block) );
                //INFO("Start = "<< mc1.genPDF(block)<<", "<<mc2.genPDF(block));
                mc1.setGenPDF(block, pdf(mc1[block], mc2[block]) / mc1.genPDF(block) );
                mc2.setGenPDF(block, pdf(mc1[block], mc2[block]) / mc2.genPDF(block) );
                }
                //maxProb = maxProb == 0 ? 1.5 * getMax(mc1) : maxProb; 
                maxProb = getMax(mc1) * 1.5;
                INFO("maxProb = "<<getMax(mc1));
                t_eval.stop();
                t_acceptReject.start(); 
                totalGenerated += mc1.size();
                //for(const auto& event1 : mc1)
                for(size_t i=0;i<mc1.size();i++)
                { 
                Event event1 = mc1[i];
                Event event2 = mc2[i];
                if ( event1.genPdf()  > maxProb || event2.genPdf() > maxProb ) {
                    std::cout << std::endl; 
                    WARNING( "PDF value exceeds norm value: " << event1.genPdf() << " > " << maxProb );
                }
                if ( event1.genPdf() > maxProb * rndm.Rndm()  ){
                    list1.push_back(event1);
                    list2.push_back(event2);
                    list1.rbegin()->setGenPdf( pdf(event1, event2) );
                    efficiencyReport[event1.index()] = true; 
                }
                else efficiencyReport[event1.index()] = false; 
                if ( list1.size() - size0 == N ) break; 
                }
                t_acceptReject.stop(); 

                // m_gps.provideEfficiencyReport( efficiencyReport );
                double efficiency = 100. * ( list1.size() - previousSize ) / (double)generatorBlock;
                pb.print( double(list1.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
                if ( list1.size() == previousSize ) {
                ERROR( "No events generated, PDF is likely to be malformed" );
                break;
                }
            } 
            pb.finish();
            t_total.stop();
            INFO("Generated " << N << " events in " << t_total << " ms");
            INFO("Generating phase space : " << t_phsp         << " ms"); 
            INFO("Evaluating PDF         : " << t_eval         << " ms"); 
            INFO("Accept/reject          : " << t_acceptReject << " ms"); 
            INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");

    


        }
        template <typename eventList_t, typename pdf_t> void generateCKM( pdf_t& pdf, eventList_t& list, const size_t& N, size_t seed, size_t generatorBlock )
        {    
            TRandom3 rndm(seed);
            //  double maxProb   = m_normalise ? 0 : 1;
            double maxProb = 1;
            auto size0       = list.size();
            double totalGenerated = 0; 
            Generator<> g(list.eventType());
        
            //pdf.reset( true );
            ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
            ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
            std::vector<bool> efficiencyReport(generatorBlock,false); 
            EventList mc( list.eventType() );
            INFO("Start generation of "<<N<<" "<<mc.eventType());
            while ( list.size() - size0 < N ) {
                t_phsp.start();
                g.fillEventListPhaseSpace(mc, generatorBlock);
                t_phsp.stop();
                t_eval.start();
            //    pdf.setEvents( mc1 );
            //    pdf.prepare();
                auto previousSize = list.size();

                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for ( size_t block=0; block < mc.nBlocks(); ++block )
                { 
                //mc1.setGenPDF(block, pdf(mc1.block(block), block) / mc1.genPDF(block) );
                //INFO("Start = "<< mc1.genPDF(block)<<", "<<mc2.genPDF(block));
                mc.setGenPDF(block, pdf(mc[block]) / mc.genPDF(block) );

                }
                //maxProb = maxProb == 0 ? 1.5 * getMax(mc1) : maxProb; 
                maxProb = getMax(mc) * 1.5;
                INFO("maxProb = "<<getMax(mc));
                t_eval.stop();
                t_acceptReject.start(); 
                totalGenerated += mc.size();
                //for(const auto& event1 : mc1)
                for(size_t i=0;i<mc.size();i++)
                { 
                Event event = mc[i];

                if ( event.genPdf()  > maxProb) {
                    std::cout << std::endl; 
                    WARNING( "PDF value exceeds norm value: " << event.genPdf() << " > " << maxProb );
                }
                if ( event.genPdf() > maxProb * rndm.Rndm()  ){
                    list.push_back(event);
                    
                    list.rbegin()->setGenPdf( pdf(event) );
                    efficiencyReport[event.index()] = true; 
                }
                else efficiencyReport[event.index()] = false; 
                if ( list.size() - size0 == N ) break; 
                }
                t_acceptReject.stop(); 

                // m_gps.provideEfficiencyReport( efficiencyReport );
                double efficiency = 100. * ( list.size() - previousSize ) / (double)generatorBlock;
                pb.print( double(list.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
                if ( list.size() == previousSize ) {
                ERROR( "No events generated, PDF is likely to be malformed" );
                break;
                }
            } 
            pb.finish();
            t_total.stop();
            INFO("Generated " << N << " events in " << t_total << " ms");
            INFO("Generating phase space : " << t_phsp         << " ms"); 
            INFO("Evaluating PDF         : " << t_eval         << " ms"); 
            INFO("Accept/reject          : " << t_acceptReject << " ms"); 
            INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");
        }
        void boostPsi3770Events(const EventList_type& list1, const EventList_type& list2, EventList& new_list1, EventList& new_list2){
//            if (boostEvents){

            TRandom3 rndm(0);
            real_t mother_mass(NamedParameter<real_t>("boostTwoBody::motherMass", 3.773));
            real_t daughter1_mass(NamedParameter<real_t>("boostTwoBody::daughter1Mass", 1.864));
            real_t daughter2_mass(NamedParameter<real_t>("boostTwoBody::daughter2Mass", 1.864));
            real_t E1 = (std::pow(mother_mass, 2) + std::pow(daughter1_mass, 2) - std::pow(daughter2_mass, 2))/(2 * mother_mass);
            real_t E2 = (std::pow(mother_mass, 2) + std::pow(daughter2_mass, 2) - std::pow(daughter1_mass, 2))/(2 * mother_mass);
            real_t theta12(rndm.Rndm() * M_PI);
            real_t phi12(rndm.Rndm() * 2 * M_PI);
            real_t gamma1(E1/daughter1_mass);
            real_t gamma2(E2/daughter2_mass);
            real_t v1(std::sqrt(1 - 1/(std::pow(gamma1, 2))));
            real_t v2(std::sqrt(1 - 1/(std::pow(gamma2, 2))));
            real_t nx1(std::sin(theta12) * std::cos(phi12));
            real_t nx2(-nx1);
            real_t ny1(std::sin(theta12) * std::sin(phi12));
            real_t ny2(-ny1);
            real_t nz1(std::cos(phi12));
            real_t nz2(-nz2);
            std::tuple<real_t, real_t, real_t> n1({nx1, ny1, nz1});
            std::tuple<real_t, real_t, real_t> n2({nx2, ny2, nz2});

            for (unsigned i=0;i<list1.size();++i){
                Event evt1(list1[i]);
                Event evt2(list2[i]);
                boost(evt1, n1, v1);
                boost(evt2, n2, v2);
                //list1[i] = evt1;
                //list2[i] = evt2;
                new_list1.push_back(evt1);
                new_list2.push_back(evt2);
            }




        }


        void boostCKMEvents(const EventList_type& list, EventList& new_list){
//            if (boostEvents){

            TRandom3 rndm(0);
            real_t mother_mass(NamedParameter<real_t>("boostTwoBody::motherMass", 5.27934));
            real_t daughter1_mass(NamedParameter<real_t>("boostTwoBody::daughter1Mass", 1.864));
            real_t daughter2_mass(NamedParameter<real_t>("boostTwoBody::daughter2Mass", 0.493677));
            real_t E1 = (std::pow(mother_mass, 2) + std::pow(daughter1_mass, 2) - std::pow(daughter2_mass, 2))/(2 * mother_mass);
            real_t E2 = (std::pow(mother_mass, 2) + std::pow(daughter2_mass, 2) - std::pow(daughter1_mass, 2))/(2 * mother_mass);
            real_t theta12(rndm.Rndm() * M_PI);
            real_t phi12(rndm.Rndm() * 2 * M_PI);
            real_t gamma1(E1/daughter1_mass);
            real_t gamma2(E2/daughter2_mass);
            real_t v1(std::sqrt(1 - 1/(std::pow(gamma1, 2))));
            real_t v2(std::sqrt(1 - 1/(std::pow(gamma2, 2))));
            real_t nx1(std::sin(theta12) * std::cos(phi12));
            real_t nx2(-nx1);
            real_t ny1(std::sin(theta12) * std::sin(phi12));
            real_t ny2(-ny1);
            real_t nz1(std::cos(phi12));
            real_t nz2(-nz2);
            std::tuple<real_t, real_t, real_t> n1({nx1, ny1, nz1});
            std::tuple<real_t, real_t, real_t> n2({nx2, ny2, nz2});



            for (unsigned i=0;i<list.size();++i){
                Event evt(list[i]);
                boost(evt, n1, v1);
                //boost(list2[i], n2, v2);
                //list[i] = evt;
                new_list.push_back(evt);
            }

        }





        void writeEvents(EventList& list, const std::string& prefix, size_t nBins=100){
            //f->cd();
            list.tree(prefix.c_str())->Write();
            auto plots = list.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack)); 
            for ( auto& plot : plots ) plot->Write((prefix + "_" + plot->GetName()).c_str());
            if( NamedParameter<bool>("plots_2d",true) == true ){
                auto proj = list.eventType().defaultProjections(nBins);
                for( size_t i = 0 ; i < proj.size(); ++i ){
                    for( size_t j = i+1 ; j < proj.size(); ++j ){ 
                        list.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack), PlotOptions::Prefix((prefix + "_" + proj[i].name() + "_" + proj[j].name()).c_str())  )->Write( ); 
                    }
                }
            }


        } 

        //template <typename amp_t=std::function<real_t(Event&)> > void writeValues(EventList_type& list, amp_t& psi, const std::string& name){
        template <typename amp_t, typename events_t> void writeValues(events_t& list, amp_t& psi, const std::string& name){
            std::string treeName =  name;
            TTree  t(treeName.c_str(), treeName.c_str());

            double x;
            auto b(t.Branch("x", &x));
            for (auto& evt : list){
                x = psi(evt);
                t.Fill();
            }
            t.Write();

        }

        void writeDalitz(EventList_type& list){
            auto s01 = [](Event& evt){ return evt.s(0, 1);};
            auto s02 = [](Event& evt){ return evt.s(0, 2);};
            auto s12 = [](Event& evt){ return evt.s(1, 2);};
            auto theta = [](Event& evt, size_t i, size_t j){ 
                double px_i = evt[4*i+0];
                double py_i = evt[4*i+1];
                double pz_i = evt[4*i+2];
                double px_j = evt[4*j+0];
                double py_j = evt[4*j+1];
                double pz_j = evt[4*j+2];
                double dot_ij = px_i * px_j + py_i * py_j + pz_i * pz_j;
                double p_i = std::pow(px_i * px_i + py_i * py_i + pz_i * pz_i, 0.5);
                double p_j = std::pow(px_j * px_j + py_j * py_j + pz_j * pz_j, 0.5);
                return std::acos(dot_ij/(p_i * p_j));
            } ;
            TTree t("Dalitz", "Dalitz");
            double x, y, z, theta12, theta13, theta23;
            auto bs01(t.Branch("s01", &x));
            auto bs02(t.Branch("s02", &y));
            auto bs12(t.Branch("s12", &z));
            auto b_theta12(t.Branch("theta_12", &theta12));
            auto b_theta13(t.Branch("theta_13", &theta13));
            auto b_theta23(t.Branch("theta_23", &theta23));
            for (auto &evt : list){
                x=s01(evt);
                y=s02(evt);
                z=s12(evt);
                theta12 = theta(evt, 1, 2);
                theta13 = theta(evt, 1, 3);
                theta23 = theta(evt, 2, 3);
                t.Fill();
            }
            t.Write();
            
        }
        

        //real_t sinTerm(std::map<std::string, std::vector<real_t> >& amps, MinuitParameterSet& MPS, std::map<std::string, ce>& cPij, EventList_type& mcSig)
 
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        std::function<real_t()> LLPsi3770(EventList& datSig, EventList& datTag, MinuitParameterSet& MPS, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, bool sameType,
                                          std::map<std::string, ce>& cPij, EventList_type& mcSig, EventList_type& mcTag){
            
            
            std::map<std::string, std::vector<real_t> > ampsDatSig(AmpArrays(datSig, A, Abar));
            std::map<std::string, std::vector<real_t> > ampsDatTag(AmpArrays(datTag, B, Bbar));
            std::map<std::string, std::vector<real_t> > ampsMCSig(AmpArrays(mcSig, A, Abar));
            std::map<std::string, std::vector<real_t> > ampsMCTag(AmpArrays(mcTag, B, Bbar));
                                            

            real_t CTag(cosTermNoCorrection(ampsMCTag, mcTag));
            real_t STag(sinTermNoCorrection(ampsMCTag, mcTag));
            real_t CTestFunc(cosTerm(ampsMCSig, MPS, cPij, mcSig));
            real_t STestFunc(sinTerm(ampsMCSig, MPS, cPij, mcSig));

            INFO("Got "<<CTestFunc<<", "<<STestFunc<<" cos/sin terms");
         //        real_t probCorr_unnorm_nonInt(const Event& evt1, const Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, bool sameTag)
            
            return [&datSig, &datTag, &MPS, &ampsMCSig, &ampsMCTag, &cPij, &mcSig, &CTag, &STag, &sameType, &ampsDatSig, &ampsDatTag, &A, &Abar, &B, &Bbar]{
                INFO("Start of function");
                INFO("stored MC has "<<ampsMCSig["A"].size()<<" |A|'s");
                INFO("stored Data has "<<ampsDatSig["A"].size()<<" |A|'s");
                real_t correction1 = phaseCorrection(datSig[0], cPij, MPS);
                //INFO("Calc Cos/Sin for norm");
                //real_t C(cosTerm(&ampsMCSig, MPS, cPij, mcSig));
                

                //real_t S(sinTerm(ampsMCSig, MPS, cPij, mcSig));
                return 1;
                /*
                INFO("Calculating normalisation");
                real_t norm = 0;
                if (sameType){
                    norm = A.norm() * Bbar.norm() + Abar.norm() * B.norm() - 2 * C * CTag - 2 * S * STag;
                }
                else{
                    norm = 2*A.norm() * Abar.norm() - 2* C * C - 2 * S * S;
                }
                INFO("Calculating LL");
                real_t ll=0;
                INFO("Start of OMP Loop");
                #pragma omp parallel for reduction (+:ll)
                for (unsigned i=0;i<datSig.size();++i){
                    real_t probInt = probCorr_unnorm_Int(datSig, datTag, i, ampsDatSig, ampsDatTag, MPS, cPij, sameType);
                    real_t probNonInt = probCorr_unnorm_nonInt(i, ampsDatSig, ampsDatTag, sameType);
                    real_t prob = probInt + probNonInt;
                    ll += log(prob/norm);
                }
                
                return -2 * ll;
                */
            } ;
            
        }
        template <class ce=CompiledExpression<real_t(const real_t*, const real_t*)> >
        std::function<real_t()> LLCKM(EventList& dat, MinuitParameterSet& MPS, CoherentSum& A, CoherentSum& Abar, int gammaSign,
                                          std::map<std::string, ce>& cPij, EventList_type& mc){
            
            std::map<std::string, std::vector<real_t> > ampsDat(AmpArrays(dat, A, Abar));
            std::map<std::string, std::vector<real_t> > ampsMC(AmpArrays(mc, A, Abar));
            return [](){
                return 1;
            };
            /*
            return [&dat, &MPS, &ampsMC, &cPij, &mc, &gammaSign, &ampsDat, &A, &Abar]{

                complex_t zB(ckm_zB(MPS, gammaSign));
                //real_t C(cosTerm(ampsMC, MPS, cPij, mc));
                //real_t S(sinTerm(ampsMC, MPS, cPij, mc));
                real_t norm(ckm_norm(A, Abar, MPS, mc, cPij, gammaSign, ampsMC));
                real_t ll=0;
                #pragma omp parallel for reduction (+:ll)
                for (unsigned i=0;i<dat.size();++i){
                    real_t probInt = probCKM_unnorm_Int(dat, i, ampsDat, MPS, cPij, gammaSign);
                    real_t probNonInt = probCKM_unnorm_nonInt(i, ampsDat, MPS, gammaSign);
                    real_t prob = probInt + probNonInt;
                    ll += log(prob/norm);
                }
                return -2 * ll;
            } ;
            */
        }


        //chi2 for correlated
        //
        struct Moment {
          double x;
          double xx;
          double N;
          std::vector<double> values;
          Moment() : x( 0 ), xx( 0 ), N( 0 ) {}
          void add( const double& value )
          {
            x += value;
            xx += value * value;
            N++;
            values.push_back( value );
          }
          void rescale( const double& val )
          {
            x *= val;
            xx *= ( val * val );
          }
          double val() { return x; }
          double var() { return N == 0 ? 0 : xx; }
        };


        
        
        template <typename pdf_t> void do_chi2_corr(EventList& datSig, EventList& datTag, EventList_type& mcSig, EventList_type& mcTag, pdf_t& pdf, MinuitParameterSet& MPS, real_t& chi2, size_t& nBins){
//            real_t chi2=0;
            size_t dofSig(datSig.eventType().dof());
//            size_t dofTag(datTag.eventType().dof());
//            INFO("dofTag = "<<dofTag);
//            if (dofTag<=0 || dofTag > 100) dofTag = 0;
//            size_t minEvents = datSig.size()/100;
//            if (minEvents < 15) minEvents = 15;
            BinDT binnedDataSig(datSig, MinEvents(NamedParameter<size_t>("QMI::corrchi2::minEvents", 15)), Dim(dofSig));
//            BinDT binnedDataSig(datSig, Dim(dofSig));
//            BinDT binnedDataTag(datTag, MinEvents(NamedParameter<size_t>("minEvents", 15)));//, Dim(dofTag));
//
//            std::map<std::pair<size_t, size_t>, real_t> dataMean; 
//            std::map<std::pair<size_t, size_t>, real_t> dataVar; 
//            std::map<std::pair<size_t, size_t>, real_t> mcMean; 
//            std::map<std::pair<size_t, size_t>, real_t> mcVar; 
//
            std::vector<Moment> data(binnedDataSig.size()); 

            std::vector<Moment> mc(binnedDataSig.size()); 
            std::vector<size_t> nDataPerBin(binnedDataSig.size());
            std::vector<size_t> nMCPerBin(binnedDataSig.size());


     




            real_t totalDatWeight(0);            
            real_t totalIntWeight(0);            
//            for (unsigned i=0;i<binnedDataSig.size();++i){
//                for (unsigned j=0;j<binnedDataTag.size();++j){
//                    std::pair<size_t, size_t> binPair({i, j});
//                    std::pair<std::pair<size_t, size_t>, real_t> entry({binPair,0});
//                    dataMean.insert(entry);
//                    dataVar.insert(entry);
//                    mcMean.insert(entry);
//                    mcVar.insert(entry);
//                }
//            }
//
            INFO("Getting data entries for "<<datSig.size()<<" events");

            real_t maxpDat =0;
            for (unsigned i=0;i<datSig.size();++i){ 
                size_t binNumSig(binnedDataSig.getBinNumber(datSig[i]));
//                INFO("At "<<i<<" out of "<<datSig.size()<<" binnum = "<<binNumSig<<" out of "<<dataMean.size());
//                size_t binNumTag(binnedDataTag.getBinNumber(datTag[i]));
//                std::pair<size_t, size_t> binPair({binNumSig, binNumTag});
//                real_t p = pdf(datSig[i], datTag[i]);///(real_t)datSig.size(); //INFO("p "<<i<<" = "<<p);
                real_t p = 1;//pdf(datSig[i], datTag[i]);///(real_t)datSig.size(); //INFO("p "<<i<<" = "<<p);
//                INFO("p = "<<p);
//                real_t p = 1;
                data[binNumSig].add(p);
                nDataPerBin[binNumSig]++;
//                dataMean[binPair] += p; 
//                dataVar[binNumSig] += p * p; 
//                dataVar[binPair] += p * p; 
               if (data[binNumSig].val()>maxpDat) maxpDat = data[binNumSig].val();
                if (i%1000 == 0){            
                auto evtWeightSig = datSig[i].weight();
                auto evtWeightTag = datTag[i].weight();
                auto evtGenPDFSig = datSig[i].genPdf();
                auto evtGenPDFTag = datTag[i].genPdf();
           //     INFO("evt ("<<i<<") = ("<<datSig[i].s(0,1)<<", "<<datSig[i].s(0, 2)<<") : bin, p, weight, genpdf = "<<binNumSig<<", "<<p<<", "<<evtWeightSig<<", "<<evtGenPDFSig<<"Tag: "<<evtWeightTag<<", "<<evtGenPDFTag);
                }
                
               
                totalDatWeight += p;
            } 
            //INFO("max p(data) = "<<maxpDat);

            real_t maxpMC =0;
            INFO("Getting mc entries");
            for (unsigned i=0;i<mcSig.size();++i){
                size_t binNumSig(binnedDataSig.getBinNumber(mcSig[i]));
//                size_t binNumTag(binnedDataTag.getBinNumber(mcTag[i]));
                real_t p= pdf(mcSig[i], mcTag[i]);///(real_t)mcSig.size(); //INFO("p "<<i<<" = "<<p);

                if (i%1000 == 0){            
                auto evtWeightSig = mcSig[i].weight();
                auto evtWeightTag = mcTag[i].weight();
                auto evtGenPDFSig = mcSig[i].genPdf();
                auto evtGenPDFTag = mcTag[i].genPdf();
                //INFO("evt "<<i<<"p weight, genpdf = "<<p<<", "<<evtWeightSig<<", "<<evtGenPDFSig<<"Tag: "<<evtWeightTag<<", "<<evtGenPDFTag);

//                INFO("evt ("<<i<<") = ("<<mcSig[i].s(0,1)<<", "<<mcSig[i].s(0, 2)<<") : bin, p, weight, genpdf = "<<binNumSig<<", "<<p<<", "<<evtWeightSig<<", "<<evtGenPDFSig<<"Tag: "<<evtWeightTag<<", "<<evtGenPDFTag);
                }

                //real_t p= 1;
 //               std::pair<size_t, size_t> binPair({binNumSig, binNumTag});
                mc[binNumSig].add(p);
                nMCPerBin[binNumSig]++;
//                mcMean[binPair] += p; 
//                mcVar[binNumSig] += p * p;  
//                mcVar[binPair] += p * p;  
                totalIntWeight += p;
               
              if (mc[binNumSig].val()>maxpMC) maxpMC = mc[binNumSig].val();

            }
            INFO("Total DatWeight = "<<totalDatWeight);
            INFO("Total IntWeight = "<<totalIntWeight);
            
            INFO("Average psi(data) = "<<totalDatWeight/(real_t)datSig.size());
            INFO("Average psi(mc) = "<<totalIntWeight/(real_t)mcSig.size());

            INFO("max p(MC) = "<<maxpMC);
            chi2 = 0;
            nBins =binnedDataSig.size() ;
            maxpMC = 0;
            real_t sumMC = 0;
            real_t sumData = 0;
            for (unsigned i=0;i<data.size();++i){
//            for (auto p : mcMean)
//                nBins++;
//                std::pair<size_t, size_t> binPair(i);
            //    mcMean[i] *= totalDatWeight/totalIntWeight;
            //    mcVar[i] *= std::pow(totalDatWeight/totalIntWeight, 2);
             //   data[i].rescale(nDataPerBin[i]);
             //   mc[i].rescale(nMCPerBin[i]);
                //mc[i].rescale
                //data[i].rescale(1/(real_t)datSig.size());
                mc[i].rescale(totalDatWeight/totalIntWeight);
                //mc[i].rescale((real_t)datSig.size()/(real_t)mcSig.size());
                sumMC += mc[i].val();
                sumData += data[i].val();

               

                if (mc[i].val() > maxpMC) maxpMC = mc[i].val();
                real_t delta = data[i].val() - mc[i].val();
                
            //    INFO("datam = "<<dataMean[i]);
             //   INFO("mcm = "<<mcMean[i]);

                real_t err = std::pow(data[i].val(), 2) + std::pow(mc[i].val(), 2);
                //real_t err = std::pow(mc[i].val(), 2);
                real_t err2 = data[i].var() + mc[i].var();
                //real_t err2 = mc[i].var();
                real_t err3 = std::pow(std::abs(err - err2), 0.5);


                chi2 += delta * delta/err2;
               // INFO("delta = "<<delta<<" data  = "<<data[i].val()<<" mc = "<<mc[i].val() << " err = "<<err<<" delta/err = "<<delta/err<<"delta2/err3 = "<<delta*delta/err3<<" err2 = "<<err2<<" err3 = "<<err3<<" chi2 = "<<chi2);
                
            }
            INFO("max MC after rescale = "<<maxpMC);
            INFO("chi2 = "<<chi2);
            INFO("sumMC = "<<sumMC<<", sumData = "<<sumData);


        }

        

        std::map<std::pair<std::string, std::string>, double> invCovarianceMatrix(std::string& logfile){
            
            
            std::ifstream infile(logfile);
            std::string line;
            std::vector<std::string> lines;
            if (infile.is_open()){
                while(std::getline(infile, line)){
                    lines.push_back(line);
                }
            }
            infile.close();
            std::string start("CovarianceMatrix");
            std::string end("ParamsOrder");
            int count = 0;
            int idxStart, idxEnd;
            std::vector<std::string> myLines;
            for (auto s : lines){
                INFO("s = "<<s);
                std::vector<std::string> candLines;
                boost::algorithm::split(candLines, s, boost::is_any_of("\t, "), boost::token_compress_on);
                bool atStart = std::find(candLines.begin(), candLines.end(), start) != candLines.end();
                bool atEnd = std::find(candLines.begin(), candLines.end(), end) != candLines.end();
                if (atStart) idxStart = count;
                if (atEnd) idxEnd = count;
                count++;
            }

            std::map<std::pair<std::string, std::string>, double> m;
            std::vector<std::string> paramNames;
            boost::algorithm::split(paramNames, lines[idxEnd], boost::is_any_of("\t, "), boost::token_compress_on);
            std::vector<std::string> paramNames2;
            INFO("paramNames has "<<paramNames.size());
            for (auto p : paramNames){
                if (p != "Free" && p!="Fix" && p != "ParamsOrder" && p != "") {
                    paramNames2.push_back(p);
                    INFO("param = "<<p);
                }
            }
            TMatrixD matrix(paramNames2.size() , paramNames2.size());
//            int idxMatrix_i = 0;
            for (int i=0;i<paramNames2.size();++i){
                int row = idxStart + i + 1;
                std::vector<std::string> rowVect;
                boost::algorithm::split(rowVect, lines[row], boost::is_any_of("\t, "), boost::token_compress_on);
                std::string param_i = paramNames2[i];
                std::cout<<"i = "<<i<<" p["<<i<<"] = "<<param_i<<"\n";
 //               int idxMatrix_j = 0;
                for (int j=0;j<paramNames2.size();++j){
                    std::string param_j = paramNames2[j];
                    int column = j ;
                    std::string param = rowVect[column];
                    //std::cout<<param_i<<" "<<param_j<<" = "<<param<<"\n";

                    double cov_param = std::stod(param);
                    matrix[i][j] = cov_param;
                }
            }
            matrix.Print();

            TMatrixD inv_matrix(paramNames2.size(), paramNames2.size());
            if (paramNames2.size() == 1){
                inv_matrix[0][0] = 1/matrix[0][0];
            }
            else{
                TDecompSVD svd(matrix);
                TMatrixD svdI = svd.Invert();
                for (int i=0;i<paramNames2.size();i++){
                    for (int j=0;j<paramNames2.size();j++){
                        inv_matrix[i][j] = svdI[i][j];
                    }

                }

//                TMatrixD inv_matrix = svd.Invert();

//                inv_matrix = matrix.Invert();

            }
            inv_matrix.Print();
            for (int i=0;i<paramNames2.size();++i){
                for (int j=0;j<paramNames2.size();++j){
                    INFO("pair = ("<<paramNames2[i]<<", "<<paramNames2[j]<<") = "<<inv_matrix[i][j]);
                    std::pair<std::string, std::string> paramPair({paramNames2[i], paramNames2[j]});
                    INFO("Making entry");
                    std::pair<std::pair<std::string, std::string>, double> cov_entry({paramPair, inv_matrix[i][j]});
                    INFO("Insert Entry");
                    m.insert(cov_entry);

                }
            }
            INFO("Done invCov");
            return m;
        }

        std::map<std::string, std::pair<double,double> > fitValAndErr(std::string& logFile){
            std::ifstream infile(logFile);
            std::string line;
            std::vector<std::string> lines;
            bool debug = NamedParameter<bool>("QMI::BDKFit::ReadPsi3770", false);
            
            if (infile.is_open()){
                while(std::getline(infile, line)){
                    lines.push_back(line);
                }
            }
            std::map<std::string, std::pair<double, double> > my_fitValAndErr;
            for (auto s : lines){
                std::vector<std::string> lineVect;
                std::vector<std::string> lineVect2;
                boost::algorithm::split(lineVect, s, boost::is_any_of("\t, "), boost::token_compress_on);
//                bool haveParameter = std::find(lineVect.start(), lineVect.end(), "Parameter") != lineVect.end();
                bool haveParameter = lineVect[0] == "Parameter";
                if (haveParameter){
                    if (debug){
                        INFO("s = "<<s);
                        INFO("len(a) = "<<lineVect.size());
                        INFO(lineVect[0]<<" "<<lineVect[1]<<" "<<lineVect[2]<<" "<<lineVect[3]<<" "<<lineVect[4]<<" "<<lineVect[5]);
                    }
                    std::string paramName = lineVect[1];
                    
                    std::string flag = lineVect[2];
                    double mean = std::stod(lineVect[3]);
                    double err = std::stod(lineVect[4]);
                    std::pair<double, double> meanAndErr({mean, err});
                    std::pair<std::string, std::pair<double, double> > my_entry({paramName, meanAndErr});
                    my_fitValAndErr.insert(my_entry);

                }

            }
            return my_fitValAndErr;
            
        }


        template<typename meanErrDict, typename invCovDict> real_t myGaussConstraint(meanErrDict meanAndErr, invCovDict iCovDict, MinuitParameterSet& MPS){
            //std::map<std::string, double> zTilde;
            //INFO("At Gauss Constraint");
            real_t my_constraint = 0;
            for (auto p1 : meanAndErr){
                //INFO("name = "<<p1.first);
                double my_z = MPS[p1.first]->mean() - meanAndErr[p1.first].first;
                double my_zTilde=0;
                for (auto p2 : meanAndErr){
                    std::pair<std::string, std::string> pairParam({p1.first, p2.first});
                    my_zTilde += (MPS[p2.first]->mean() - meanAndErr[p2.first].first) * iCovDict[pairParam];
                }
                my_constraint += my_z * my_zTilde;
            
            }
            return my_constraint;
        }


        namespace QMIPlotOptions{
            DECLARE_ARGUMENT(minX, double);
            DECLARE_ARGUMENT(maxX, double);
            DECLARE_ARGUMENT(minY, double);
            DECLARE_ARGUMENT(maxY, double);
            DECLARE_ARGUMENT(name, std::string);
            DECLARE_ARGUMENT(xAxisTitle, std::string);
            DECLARE_ARGUMENT(yAxisTitle, std::string);
            DECLARE_ARGUMENT(units, std::string);
            DECLARE_ARGUMENT(posXFunction, std::function<real_t(Event)>);
            DECLARE_ARGUMENT(posYFunction, std::function<real_t(Event)>);

            

        }

//template <typename pdf_t> std::tuple<TH1D*, THStack*> proj_1D_psi3770(const EventList_type& events_sig, const EventList_type& events_tag, const pdf_t& weightFunction, const ArgumentPack& args) const
        template <typename pdf_t> TH1D* proj_1D_psi3770(const EventList_type& events_sig, const EventList_type& events_tag, const pdf_t& weightFunction, const ArgumentPack& args) 
        {
       // TH1D* hist; 
        double norm_sum = args.getArg<PlotOptions::Norm>(1).val;
        std::string prefix = args.getArg<PlotOptions::Prefix>().val;
        bool autowrite     = args.get<PlotOptions::AutoWrite>() != nullptr;
        //  THStack* stack     = args.getArg<PlotOptions::AddTo>(new THStack()).val;
        auto selection      = args.getArg<PlotOptions::Selection>().val;
        auto my_units = args.getArg<QMIPlotOptions::units>("").val;
        auto my_nBins = args.getArg<PlotOptions::Bins>(100).val;
        auto my_min = args.getArg<QMIPlotOptions::minX>(0).val;
        auto my_max = args.getArg<QMIPlotOptions::maxX>(3).val;
        auto my_name = args.getArg<QMIPlotOptions::name>("s01").val;
        auto my_width = (my_max - my_min)/(real_t)my_nBins;
        auto my_posFunction = args.getArg<QMIPlotOptions::posXFunction>([](Event evt){return evt.s(0, 1);}).val;
        auto my_xAxisTitle = args.getArg<QMIPlotOptions::xAxisTitle>("s01").val;

        if( prefix != "" ) prefix = prefix +"_";

            std::string p = ( prefix==""?"":prefix+"_");
            TH1D* plot = new TH1D( (p + my_name ).c_str(),"", my_nBins, my_min, my_max);
            if( my_units != "" ){ 
                plot->GetXaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", my_xAxisTitle.c_str(),my_units.c_str() ).c_str() );
                plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f %s)", my_width,my_units.c_str()).c_str() );
            }
            else {
                plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f)",my_width).c_str() );
                plot->GetXaxis()->SetTitle( my_xAxisTitle.c_str() );
            }
            plot->GetYaxis()->SetTitleOffset(1.35);
            plot->SetMarkerSize(0);
            plot->SetMinimum(0);
        
            INFO("Returning plot: [" << my_min << " " << my_max << "] " << my_name << " " << 
            plot->GetXaxis()->GetBinLowEdge(1) << " " <<
            plot->GetXaxis()->GetBinLowEdge(1 + my_nBins)
            );




        //for( const auto& evt : events )
        for( unsigned i=0;i<events_sig.size();++i ){
            Event evt = events_sig[i];
            if( selection != nullptr && !selection(evt) ) continue;
        //    auto pos = operator()(evt);
            auto pos = my_posFunction(evt);
            auto weight = weightFunction(events_sig[i], events_tag[i]);
            plot->Fill( pos, evt.weight() * weight / evt.genPdf() ); 
        }
        //std::sort( std::begin(plots), std::end(plots), [](auto& h1, auto& h2){ return h1->Integral() < h2->Integral() ; } );
        //double total = std::accumulate( std::begin(plots), std::end(plots), 0.0, [](double& t, auto& h){ return t + h->Integral() ; } ); 
        double total = plot->Integral();
        
        if( norm_sum != -1 )
        {
            if( total == 0 ) ERROR("Norm = " << total );
        //    else for( auto& h : plots ) h->Scale( norm_sum / total );
            else{
                plot->Scale(norm_sum/total);
            }
        }
        //stack->SetName( (prefix + name() + "_stack").c_str());
        //  for( auto& h : hists ){
        //   stack->Add(h, "C HIST");
        //   if( autowrite ) h->Write();
        // }
        // stack->Add(hist, "C_HIST");
        
        //if( autowrite ) {
        //    hist->Write();
        //    stack->Write();
        //}
        //return {hist, stack};
        return plot;
        }

        template <typename pdf_t> TH2D* proj_2D_psi3770(const EventList_type& events_sig, const EventList_type& events_tag, const pdf_t& weightFunction, const ArgumentPack& args){
        //TH2D* hist; 
        double norm_sum = args.getArg<PlotOptions::Norm>(1).val;
        std::string prefix = args.getArg<PlotOptions::Prefix>().val;
        bool autowrite     = args.get<PlotOptions::AutoWrite>() != nullptr;
        //  THStack* stack     = args.getArg<PlotOptions::AddTo>(new THStack()).val;
        auto selection      = args.getArg<PlotOptions::Selection>().val;
        auto my_units = args.getArg<QMIPlotOptions::units>("").val;
        auto my_nBins = args.getArg<PlotOptions::Bins>(100).val;
        auto my_minX = args.getArg<QMIPlotOptions::minX>(0).val;
        auto my_maxX = args.getArg<QMIPlotOptions::maxX>(3).val;
        auto my_minY = args.getArg<QMIPlotOptions::minY>(0).val;
        auto my_maxY = args.getArg<QMIPlotOptions::maxY>(3).val;

        auto my_name = args.getArg<QMIPlotOptions::name>("s01_vs_s02").val;
        auto my_width = (my_maxX - my_minX)/(real_t)my_nBins;
        auto my_posXFunction = args.getArg<QMIPlotOptions::posXFunction>([](Event evt){return evt.s(0, 1);}).val;
        auto my_posYFunction = args.getArg<QMIPlotOptions::posYFunction>([](Event evt){return evt.s(0, 2);}).val;
        auto my_xAxisTitle = args.getArg<QMIPlotOptions::xAxisTitle>("s01").val;
        auto my_yAxisTitle = args.getArg<QMIPlotOptions::yAxisTitle>("s02").val;

        if( prefix != "" ) prefix = prefix +"_";
        
        std::string p = ( prefix==""?"":prefix+"_");
        TH2D* plot = new TH2D( (p + my_name ).c_str(),"", (size_t)std::pow(my_nBins, 0.5), my_minX, my_maxX, (size_t)std::pow(my_nBins, 0.5), my_maxY, my_minY);
        if( my_units != "" ){ 
            plot->GetXaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", my_xAxisTitle.c_str(),my_units.c_str() ).c_str() );
            plot->GetYaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", my_yAxisTitle.c_str(),my_units.c_str() ).c_str() );
    //        plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f %s)", my_width,units.c_str()).c_str() );
        }
        else {
            //plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f)",my_width).c_str() );
            plot->GetXaxis()->SetTitle( my_xAxisTitle.c_str() );
            plot->GetYaxis()->SetTitle( my_yAxisTitle.c_str() );
        }
        plot->GetYaxis()->SetTitleOffset(1.35);
        plot->SetMarkerSize(0);
        plot->SetMinimum(0);
    
        INFO("Returning plot: [" << my_minX << " " << my_maxX << "] " << my_name << " " << 
        plot->GetXaxis()->GetBinLowEdge(1) << " " <<
        plot->GetXaxis()->GetBinLowEdge(1 + my_nBins)
        );

        //for( const auto& evt : events ){
        for( unsigned i=0;i<events_sig.size();++i ){
            Event evt = events_sig[i];
            if( selection != nullptr && !selection(evt) ) continue;
        //    auto pos = operator()(evt);
            auto posX = my_posXFunction(evt);
            auto posY = my_posYFunction(evt);
            auto weight = weightFunction(events_sig[i], events_tag[i]);
            plot->Fill( posX, posY, evt.weight() * weight / evt.genPdf() ); 
        }
        //std::sort( std::begin(hists), std::end(hists), [](auto& h1, auto& h2){ return h1->Integral() < h2->Integral() ; } );
        //double total = std::accumulate( std::begin(hists), std::end(hists), 0.0, [](double& t, auto& h){ return t + h->Integral() ; } ); 
        double total = plot->Integral();
        
        if( norm_sum != -1 )
        {
            if( total == 0 ) ERROR("Norm = " << total );
        //    else for( auto& h : hists ) h->Scale( norm_sum / total );
            else{
                plot->Scale(norm_sum/total);
            }
        } 
        
        return plot;
        }

        template <typename pdf_t> TH2D* proj_2D_BDK(const EventList_type& events, const pdf_t& weightFunction, const ArgumentPack& args){
        //TH2D* hist; 
        double norm_sum = args.getArg<PlotOptions::Norm>(1).val;
        std::string prefix = args.getArg<PlotOptions::Prefix>().val;
        bool autowrite     = args.get<PlotOptions::AutoWrite>() != nullptr;
        //  THStack* stack     = args.getArg<PlotOptions::AddTo>(new THStack()).val;
        auto selection      = args.getArg<PlotOptions::Selection>().val;
        auto my_units = args.getArg<QMIPlotOptions::units>("").val;
        auto my_nBins = args.getArg<PlotOptions::Bins>(100).val;
        auto my_minX = args.getArg<QMIPlotOptions::minX>(0).val;
        auto my_maxX = args.getArg<QMIPlotOptions::maxX>(3).val;
        auto my_minY = args.getArg<QMIPlotOptions::minY>(0).val;
        auto my_maxY = args.getArg<QMIPlotOptions::maxY>(3).val;

        auto my_name = args.getArg<QMIPlotOptions::name>("s01_vs_s02").val;
        auto my_width = (my_maxX - my_minX)/(real_t)my_nBins;
        auto my_posXFunction = args.getArg<QMIPlotOptions::posXFunction>([](Event evt){return evt.s(0, 1);}).val;
        auto my_posYFunction = args.getArg<QMIPlotOptions::posYFunction>([](Event evt){return evt.s(0, 2);}).val;
        auto my_xAxisTitle = args.getArg<QMIPlotOptions::xAxisTitle>("s01").val;
        auto my_yAxisTitle = args.getArg<QMIPlotOptions::yAxisTitle>("s02").val;

        if( prefix != "" ) prefix = prefix +"_";
        
        std::string p = ( prefix==""?"":prefix+"_");
        TH2D* plot = new TH2D( (p + my_name ).c_str(),"", (size_t)std::pow(my_nBins, 0.5), my_minX, my_maxX, (size_t)std::pow(my_nBins, 0.5), my_maxY, my_minY);
        if( my_units != "" ){ 
            plot->GetXaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", my_xAxisTitle.c_str(),my_units.c_str() ).c_str() );
            plot->GetYaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", my_yAxisTitle.c_str(),my_units.c_str() ).c_str() );
    //        plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f %s)", my_width,units.c_str()).c_str() );
        }
        else {
            //plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f)",my_width).c_str() );
            plot->GetXaxis()->SetTitle( my_xAxisTitle.c_str() );
            plot->GetYaxis()->SetTitle( my_yAxisTitle.c_str() );
        }
        plot->GetYaxis()->SetTitleOffset(1.35);
        plot->SetMarkerSize(0);
        plot->SetMinimum(0);
    
        INFO("Returning plot: [" << my_minX << " " << my_maxX << "] " << my_name << " " << 
        plot->GetXaxis()->GetBinLowEdge(1) << " " <<
        plot->GetXaxis()->GetBinLowEdge(1 + my_nBins)
        );

        //for( const auto& evt : events ){
        for( unsigned i=0;i<events.size();++i ){
            Event evt = events[i];
            if( selection != nullptr && !selection(evt) ) continue;
        //    auto pos = operator()(evt);
            auto posX = my_posXFunction(evt);
            auto posY = my_posYFunction(evt);
            auto weight = weightFunction(events[i]);
            plot->Fill( posX, posY, evt.weight() * weight / evt.genPdf() ); 
        }
        //std::sort( std::begin(hists), std::end(hists), [](auto& h1, auto& h2){ return h1->Integral() < h2->Integral() ; } );
        //double total = std::accumulate( std::begin(hists), std::end(hists), 0.0, [](double& t, auto& h){ return t + h->Integral() ; } ); 
        double total = plot->Integral();
        
        if( norm_sum != -1 )
        {
            if( total == 0 ) ERROR("Norm = " << total );
        //    else for( auto& h : hists ) h->Scale( norm_sum / total );
            else{
                plot->Scale(norm_sum/total);
            }
        } 
        
        return plot;
        }

    real_t logPoisson(real_t x, real_t m){
        real_t y = 0;

        //y = -(m) + (x) * std::log(m) - std::lgamma(x + 1);
        //y = -0.5 * std::pow(x - m, 2)/(m + x) - 0.5 * std::log(2 * M_PI * (x + m));
         if (x!= 0) y = -0.5 * std::pow(x - m, 2)/(x);// - 0.5 * std::log(2 * M_PI * (x + 1));
         return y;
}
    real_t chi2(real_t x, real_t m){
        real_t y = 0;

        y = (x -m) * (x-m)/(m); 
 
         return y;
}


    }//QMI Namespace

}//AmpGen Namespace
#endif
