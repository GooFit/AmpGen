#ifndef PHASECORRECTION_H_
#define PHASECORRECTION_H_
#include <functional>
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Event.h"
#include <TF1.h>
#include <Math/WrappedTF1.h>
//#include <Math/SpecFunMathMore.h>
//#include "AmpGen/Polynomials.h"
namespace AmpGen{
    class PhaseCorrection{
        public:
            real_t fastLegendre(real_t x, size_t n){
                switch (n)
                    {
                    case 0/* constant-expression */:
                        /* code */
                        return 1;
            

                    case 1:
                        return x;
                    case 2:
                        return 1.5*std::pow(x, 2) - 0.5;
                    case 3:
                        return 2.5*std::pow(x, 3) - 1.5*x;
                    case 4:
                        return 4.375*std::pow(x, 4) - 3.75*std::pow(x, 2) + 0.375;
                    case 5:
                        return 7.875*std::pow(x, 5) - 8.75*std::pow(x, 3) + 1.875*x;
                    case 6:
                        return 14.4375*std::pow(x, 6) - 19.6875*std::pow(x, 4) + 6.5625*std::pow(x, 2) - 0.3125;
                    case 7:
                        return 26.8125*std::pow(x, 7) - 43.3125*std::pow(x, 5) + 19.6875*std::pow(x, 3) - 2.1875*x;
                    case 8:
                        return 50.2734375*std::pow(x, 8) - 93.84375*std::pow(x, 6) + 54.140625*std::pow(x, 4) - 9.84375*std::pow(x, 2) + 0.2734375;
                    case 9:
                        return 94.9609375*std::pow(x, 9) - 201.09375*std::pow(x, 7) + 140.765625*std::pow(x, 5) - 36.09375*std::pow(x, 3) + 2.4609375*x;
                    case 10:
                        return 180.42578125*std::pow(x, 10) - 427.32421875*std::pow(x, 8) + 351.9140625*std::pow(x, 6) - 117.3046875*std::pow(x, 4) + 13.53515625*std::pow(x, 2) - 0.24609375;
                    case 11:
                        return 344.44921875*std::pow(x, 11) - 902.12890625*std::pow(x, 9) + 854.6484375*std::pow(x, 7) - 351.9140625*std::pow(x, 5) + 58.65234375*std::pow(x, 3) - 2.70703125*x;
                    case 12:
                        return 660.1943359375*std::pow(x, 12) - 1894.470703125*std::pow(x, 10) + 2029.7900390625*std::pow(x, 8) - 997.08984375*std::pow(x, 6) + 219.9462890625*std::pow(x, 4) - 17.595703125*std::pow(x, 2) + 0.2255859375;
                    case 13:
                        return 1269.6044921875*std::pow(x, 13) - 3961.166015625*std::pow(x, 11) + 4736.1767578125*std::pow(x, 9) - 2706.38671875*std::pow(x, 7) + 747.8173828125*std::pow(x, 5) - 87.978515625*std::pow(x, 3) + 2.9326171875*x;
                    case 14:
                        return 2448.52294921875*std::pow(x, 14) - 8252.42919921875*std::pow(x, 12) + 10893.2065429688*std::pow(x, 10) - 7104.26513671875*std::pow(x, 8) + 2368.08837890625*std::pow(x, 6) - 373.90869140625*std::pow(x, 4) + 21.99462890625*std::pow(x, 2) - 0.20947265625;
                    case 15:
                        return 4733.81103515625*std::pow(x, 15) - 17139.6606445313*std::pow(x, 13) + 24757.2875976563*std::pow(x, 11) - 18155.3442382812*std::pow(x, 9) + 7104.26513671875*std::pow(x, 7) - 1420.85302734375*std::pow(x, 5) + 124.63623046875*std::pow(x, 3) - 3.14208984375*x;
                    case 16:
                        return 9171.75888061523*std::pow(x, 16) - 35503.5827636719*std::pow(x, 14) + 55703.8970947266*std::pow(x, 12) - 45388.3605957031*std::pow(x, 10) + 20424.7622680664*std::pow(x, 8) - 4972.98559570313*std::pow(x, 6) + 592.022094726563*std::pow(x, 4) - 26.707763671875*std::pow(x, 2) + 0.196380615234375;
                    case 17:
                        return 17804.002532959*std::pow(x, 17) - 73374.0710449219*std::pow(x, 15) + 124262.539672852*std::pow(x, 13) - 111407.794189453*std::pow(x, 11) + 56735.4507446289*std::pow(x, 9) - 16339.8098144531*std::pow(x, 7) + 2486.49279785156*std::pow(x, 5) - 169.149169921875*std::pow(x, 3) + 3.33847045898438*x;
                    case 18:
                        return 34618.8938140869*std::pow(x, 18) - 151334.021530151*std::pow(x, 16) + 275152.766418457*std::pow(x, 14) - 269235.502624512*std::pow(x, 12) + 153185.717010498*std::pow(x, 10) - 51061.905670166*std::pow(x, 8) + 9531.55572509766*std::pow(x, 6) - 888.033142089844*std::pow(x, 4) + 31.7154693603516*std::pow(x, 2) - 0.185470581054688;
                    case 19:
                        return 67415.7405853271*std::pow(x, 19) - 311570.044326782*std::pow(x, 17) + 605336.086120605*std::pow(x, 15) - 642023.121643066*std::pow(x, 13) + 403853.253936768*std::pow(x, 11) - 153185.717010498*std::pow(x, 9) + 34041.2704467773*std::pow(x, 7) - 4084.95245361328*std::pow(x, 5) + 222.008285522461*std::pow(x, 3) - 3.52394104003906*x;
                    default:
                        return (2 * n - 1)/n * x * fastLegendre(x, n - 1) - (n - 1)/n * fastLegendre(x, n - 2);
                    }
            
            }
            real_t calcPoly(real_t x, real_t y, size_t order){
                real_t p=0;
                for (size_t i=0; i<order+1;i++){
                    size_t i1 = i;
                    size_t i2 = order - i;
                    p+=m_mps["pCorrelatedSum::C"+std::to_string(i1)+std::to_string(i2)]->mean() * fastLegendre(x, i1) * fastLegendre(y,2 * i2 + 1);
                    if (m_debug) INFO("f"<<i1<<i2<<" = "<<p);
                }
                return p;
                
            }

            


            PhaseCorrection(const MinuitParameterSet& mps) : m_mps(mps), m_order(NamedParameter<size_t>("pCorrelatedSum::Order", 2)), m_debug(NamedParameter<bool>("PhaseCorrection::debug", true)), m_start(NamedParameter<size_t>("pCorrelatedSum::Start", 0)) {}
            std::vector<real_t> getXY(Event event){
                 auto x = event.s(0,1);
                auto y = event.s(0,2);
                auto z = event.s(1,2);
                /*
                auto mp = sqrt(event.s(1,1))/2;
                auto mm = sqrt(event.s(2,2))/2;
                auto mK = sqrt(event.s(0,0))/2;
                //R__LOAD_LIBRARY(libMathMore);
                //double a = ROOT::Math::legendre(1,1);
                //double a = ROOT::Math::sin(1);
                //TF1 * f1;
                //double a = Math::sin(2);
                auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
                Expression xmin = pow(mp + mK, 2);
                Expression xmax = pow(mD - mm, 2);
                Expression x0 = (xmax + xmin)/2;
                Expression ymin = pow(mp + mK, 2);
                Expression ymax = pow(mD - mp, 2);
                Expression y0 = (ymax + ymin)/2;
                */
                real_t xmin = 0.406004;
                real_t xmax = 2.976556;
                real_t X = (2 * x - xmax - xmin)/(xmax - xmin);
                real_t Y = (2 * y - xmax - xmin)/(xmax - xmin);
                auto zp = (X+Y)/2;
                auto zm = (X-Y)/2;
                return std::vector<real_t>({zp,zm});

            }
            real_t calcCorr(Event event){
                auto XY=getXY(event);
                auto zp = XY[0];
                auto zm = XY[1];

                ProfileClock pc;
                pc.start();
                real_t loopTime=0;
                real_t coeffTime=0;
                real_t polyTime=0;
                real_t corr =0;
                for (auto i=0; i < m_order+1; i++){
                    real_t sum_i=0;
                    for (auto j=0; j<m_order+1-i; j++){

                        ProfileClock coeffClock;
                       
                        //INFO("Time to get C"<<i<<j<<" = "<<coeffClock.t_duration);
                        ProfileClock polyClock;
                        polyClock.start();
                        //sum_i += Cij*fcn::legendre(zp,i) * fcn::legendre(zm,2*j+1);
                        sum_i += m_mps["pCorrelatedSum::C"+std::to_string(i) + std::to_string(j)]->mean()*fastLegendre(zp,i) * fastLegendre(zm,2*j+1);
                        //sum_i += Cij*fcn::pow(zp,i) * fcn::pow(zm,2*j+1);
                        polyClock.stop();
                        polyTime += polyClock.t_duration;
               
                    }
                    corr += sum_i;
                }
                pc.stop();
                loopTime = polyTime;
                if (m_debug){
                    INFO("Time for Poly = "<<polyTime);
                    INFO("Times for Coeff + poly = "<<loopTime);
                    INFO("Time to evalulate event = "<<pc.t_duration);
                }
                return corr;
            }


            Expression fastCorr(const Event event){
                auto XY=getXY(event);
                auto x = XY[0];
                auto y = XY[1];
                switch(m_order){
                    case 0:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1);
                    case 1:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1);
                    case 2:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1);
                    case 3:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1);
                    case 4:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1);
                    case 5:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1);
                    case 6:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1);
                    case 7:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C25"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1);
                    case 8:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["pCorrelatedSum::C08"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C17"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C25"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C26"]->mean() * fastLegendre(x,2) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C35"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C44"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C53"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C62"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C71"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C80"]->mean() * fastLegendre(x,8) * fastLegendre(y,1);
                    case 9:
                        return m_mps["pCorrelatedSum::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["pCorrelatedSum::C08"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["pCorrelatedSum::C09"]->mean() * fastLegendre(x,0) * fastLegendre(y,19) + m_mps["pCorrelatedSum::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C17"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["pCorrelatedSum::C18"]->mean() * fastLegendre(x,1) * fastLegendre(y,17) + m_mps["pCorrelatedSum::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C25"]->mean() * fcn::pow(x,2) * fcn::pow(y,11) + m_mps["pCorrelatedSum::C26"]->mean() * fcn::pow(x,2) * fcn::pow(y,13) + m_mps["pCorrelatedSum::C27"]->mean() * fcn::pow(x,2) * fcn::pow(y,15) + m_mps["pCorrelatedSum::C30"]->mean() * fcn::pow(x,3) * fcn::pow(y,1) + m_mps["pCorrelatedSum::C31"]->mean() * fcn::pow(x,3) * fcn::pow(y,3) + m_mps["pCorrelatedSum::C32"]->mean() * fcn::pow(x,3) * fcn::pow(y,5) + m_mps["pCorrelatedSum::C33"]->mean() * fcn::pow(x,3) * fcn::pow(y,7) + m_mps["pCorrelatedSum::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C35"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C36"]->mean() * fastLegendre(x,3) * fastLegendre(y,13) + m_mps["pCorrelatedSum::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C44"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C45"]->mean() * fastLegendre(x,4) * fastLegendre(y,11) + m_mps["pCorrelatedSum::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C53"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C54"]->mean() * fastLegendre(x,5) * fastLegendre(y,9) + m_mps["pCorrelatedSum::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C62"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C63"]->mean() * fastLegendre(x,6) * fastLegendre(y,7) + m_mps["pCorrelatedSum::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C71"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C72"]->mean() * fastLegendre(x,7) * fastLegendre(y,5) + m_mps["pCorrelatedSum::C80"]->mean() * fastLegendre(x,8) * fastLegendre(y,1) + m_mps["pCorrelatedSum::C81"]->mean() * fastLegendre(x,8) * fastLegendre(y,3) + m_mps["pCorrelatedSum::C90"]->mean() * fastLegendre(x,9) * fastLegendre(y,1);
             
                    default:
                        return calcCorr(event);
                }



            }


            real_t calcCorrL(Event event){
                auto XY=getXY(event);
                auto x = XY[0];
                auto y = XY[1];
                real_t p = 0;
                for (size_t i=m_start;i<m_order+1 ; i++){
                   
                    p+=calcPoly(x,y,i);
                    if (m_debug) INFO("f"<<i<<" = "<<p);

                }
                return p;

            }
            std::function<real_t(real_t,real_t,size_t,MinuitParameterSet)> makefunc(){
                size_t order = m_order;
                std::function<real_t(real_t,real_t,int,MinuitParameterSet)> f = 
                    [](real_t x, real_t y,size_t order, MinuitParameterSet mps){
                        real_t corr=0; 
                        
                        for (auto i=0; i < order+1; i++){
                        real_t sum_i=0;
                            for (auto j=0; j<order+1-i; j++){
                            real_t Cij = mps["pCorrelatedSum::C"+std::to_string(i) + std::to_string(j)]->mean();
                            
                            //sum_i += Cij*fastLegendre(x,i) * fastLegendre(y,2*j+1);
                            sum_i += Cij*std::pow(x,i)* std::pow(y,2*j+1);
                            }
                        corr += sum_i;
                        }
                        
                    return corr;
                };
                return f;
            }
            


            

            
        private:
            MinuitParameterSet m_mps;
            size_t m_order;
            size_t m_start;
            bool m_debug;
          

    };
}
#endif