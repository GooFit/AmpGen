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
#include "TMatrixTSym.h"
#include <cmath>
//#include <Math/SpecFunMathMore.h>
//#include "AmpGen/Polynomials.h"
namespace AmpGen{
    class PhaseCorrection{
        public:
            const real_t fastLegendre(real_t x, size_t n) const {
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

            const size_t getOrder() const {
                return m_order;
            }

            const real_t Poly1D(real_t x, size_t i) const {

                if (m_PolyType=="antiSym_legendre"){
                    return fastLegendre(x, i);
                }
                else if (m_PolyType=="antiSym_simple"){
                    return std::pow(x, i); 
                }
                return 0;
                
            }




            const real_t Poly2D(real_t x, real_t y, size_t i, size_t j)const {
                if (m_PolyType=="antiSym_legendre"){
                    return fastLegendre(x, i) * fastLegendre(y, j);
                }
                else if (m_PolyType=="antiSym_simple"){
                    return std::pow(x, i) * std::pow(y, j);
                }
                else{
                    return 0;
                }
                return 0;
                
            }



            const real_t calcPoly(real_t x, real_t y, size_t order) const {
                real_t p=0;
                for (size_t i=0; i<order+1;i++){
                    size_t i1 = i;
                    size_t i2 = order - i;
		    if (m_PolyType=="antiSym_legendre"){
                    	p+=m_mps["PhaseCorrection::C"+std::to_string(i1)+"_"+std::to_string(2*i2+1)]->mean() * fastLegendre(x, i1) * fastLegendre(y,2 * i2 + 1);
                    	//p+=m_mps["PhaseCorrection::C"+std::to_string(i1)+"_"+std::to_string(2*i2+1)]->mean() * std::legendre(i1,x) * std::legendre(2 * i2 + 1, y);
		    }
		    else if (m_PolyType=="antiSym_simple"){
                    	p+=m_mps["PhaseCorrection::C"+std::to_string(i1)+"_"+std::to_string(2*i2+1)]->mean() * std::pow(x, i1) * std::pow(y,2 * i2 + 1);
		    }


                    if (m_debug) INFO("f"<<i1<<i2<<" = "<<p);
                }
                return p;
                
            }

            const real_t Legendre2D(real_t x, real_t y) const {
                
                switch (m_order)
                {

case 0:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y);
break;
case 1:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y);
break;
case 2:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y);
break;
case 3:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y);
break;
case 4:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y);
break;
case 5:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C0_11"]->mean() * std::legendre(0, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_9"]->mean() * std::legendre(1, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_7"]->mean() * std::legendre(2, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_5"]->mean() * std::legendre(3, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C4_3"]->mean() * std::legendre(4, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C5_1"]->mean() * std::legendre(5, x) * std::legendre(1, y);
break;
case 6:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C0_11"]->mean() * std::legendre(0, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C0_13"]->mean() * std::legendre(0, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_9"]->mean() * std::legendre(1, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C1_11"]->mean() * std::legendre(1, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_7"]->mean() * std::legendre(2, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C2_9"]->mean() * std::legendre(2, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_5"]->mean() * std::legendre(3, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C3_7"]->mean() * std::legendre(3, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C4_3"]->mean() * std::legendre(4, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C4_5"]->mean() * std::legendre(4, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C5_1"]->mean() * std::legendre(5, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C5_3"]->mean() * std::legendre(5, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C6_1"]->mean() * std::legendre(6, x) * std::legendre(1, y);
break;
case 7:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C0_11"]->mean() * std::legendre(0, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C0_13"]->mean() * std::legendre(0, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C0_15"]->mean() * std::legendre(0, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_9"]->mean() * std::legendre(1, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C1_11"]->mean() * std::legendre(1, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C1_13"]->mean() * std::legendre(1, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_7"]->mean() * std::legendre(2, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C2_9"]->mean() * std::legendre(2, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C2_11"]->mean() * std::legendre(2, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_5"]->mean() * std::legendre(3, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C3_7"]->mean() * std::legendre(3, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C3_9"]->mean() * std::legendre(3, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C4_3"]->mean() * std::legendre(4, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C4_5"]->mean() * std::legendre(4, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C4_7"]->mean() * std::legendre(4, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C5_1"]->mean() * std::legendre(5, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C5_3"]->mean() * std::legendre(5, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C5_5"]->mean() * std::legendre(5, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C6_1"]->mean() * std::legendre(6, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C6_3"]->mean() * std::legendre(6, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C7_1"]->mean() * std::legendre(7, x) * std::legendre(1, y);
break;
case 8:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C0_11"]->mean() * std::legendre(0, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C0_13"]->mean() * std::legendre(0, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C0_15"]->mean() * std::legendre(0, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C0_17"]->mean() * std::legendre(0, x) * std::legendre(17, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_9"]->mean() * std::legendre(1, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C1_11"]->mean() * std::legendre(1, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C1_13"]->mean() * std::legendre(1, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C1_15"]->mean() * std::legendre(1, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_7"]->mean() * std::legendre(2, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C2_9"]->mean() * std::legendre(2, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C2_11"]->mean() * std::legendre(2, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C2_13"]->mean() * std::legendre(2, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_5"]->mean() * std::legendre(3, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C3_7"]->mean() * std::legendre(3, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C3_9"]->mean() * std::legendre(3, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C3_11"]->mean() * std::legendre(3, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C4_3"]->mean() * std::legendre(4, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C4_5"]->mean() * std::legendre(4, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C4_7"]->mean() * std::legendre(4, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C4_9"]->mean() * std::legendre(4, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C5_1"]->mean() * std::legendre(5, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C5_3"]->mean() * std::legendre(5, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C5_5"]->mean() * std::legendre(5, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C5_7"]->mean() * std::legendre(5, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C6_1"]->mean() * std::legendre(6, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C6_3"]->mean() * std::legendre(6, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C6_5"]->mean() * std::legendre(6, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C7_1"]->mean() * std::legendre(7, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C7_3"]->mean() * std::legendre(7, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C8_1"]->mean() * std::legendre(8, x) * std::legendre(1, y);
break;
case 9:
    return m_mps["PhaseCorrection::C0_1"]->mean() * std::legendre(0, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C0_3"]->mean() * std::legendre(0, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C0_5"]->mean() * std::legendre(0, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C0_7"]->mean() * std::legendre(0, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C0_9"]->mean() * std::legendre(0, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C0_11"]->mean() * std::legendre(0, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C0_13"]->mean() * std::legendre(0, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C0_15"]->mean() * std::legendre(0, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C0_17"]->mean() * std::legendre(0, x) * std::legendre(17, y) + m_mps["PhaseCorrection::C0_19"]->mean() * std::legendre(0, x) * std::legendre(19, y) + m_mps["PhaseCorrection::C1_1"]->mean() * std::legendre(1, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C1_3"]->mean() * std::legendre(1, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C1_5"]->mean() * std::legendre(1, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C1_7"]->mean() * std::legendre(1, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C1_9"]->mean() * std::legendre(1, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C1_11"]->mean() * std::legendre(1, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C1_13"]->mean() * std::legendre(1, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C1_15"]->mean() * std::legendre(1, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C1_17"]->mean() * std::legendre(1, x) * std::legendre(17, y) + m_mps["PhaseCorrection::C2_1"]->mean() * std::legendre(2, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C2_3"]->mean() * std::legendre(2, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C2_5"]->mean() * std::legendre(2, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C2_7"]->mean() * std::legendre(2, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C2_9"]->mean() * std::legendre(2, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C2_11"]->mean() * std::legendre(2, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C2_13"]->mean() * std::legendre(2, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C2_15"]->mean() * std::legendre(2, x) * std::legendre(15, y) + m_mps["PhaseCorrection::C3_1"]->mean() * std::legendre(3, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C3_3"]->mean() * std::legendre(3, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C3_5"]->mean() * std::legendre(3, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C3_7"]->mean() * std::legendre(3, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C3_9"]->mean() * std::legendre(3, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C3_11"]->mean() * std::legendre(3, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C3_13"]->mean() * std::legendre(3, x) * std::legendre(13, y) + m_mps["PhaseCorrection::C4_1"]->mean() * std::legendre(4, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C4_3"]->mean() * std::legendre(4, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C4_5"]->mean() * std::legendre(4, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C4_7"]->mean() * std::legendre(4, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C4_9"]->mean() * std::legendre(4, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C4_11"]->mean() * std::legendre(4, x) * std::legendre(11, y) + m_mps["PhaseCorrection::C5_1"]->mean() * std::legendre(5, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C5_3"]->mean() * std::legendre(5, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C5_5"]->mean() * std::legendre(5, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C5_7"]->mean() * std::legendre(5, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C5_9"]->mean() * std::legendre(5, x) * std::legendre(9, y) + m_mps["PhaseCorrection::C6_1"]->mean() * std::legendre(6, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C6_3"]->mean() * std::legendre(6, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C6_5"]->mean() * std::legendre(6, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C6_7"]->mean() * std::legendre(6, x) * std::legendre(7, y) + m_mps["PhaseCorrection::C7_1"]->mean() * std::legendre(7, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C7_3"]->mean() * std::legendre(7, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C7_5"]->mean() * std::legendre(7, x) * std::legendre(5, y) + m_mps["PhaseCorrection::C8_1"]->mean() * std::legendre(8, x) * std::legendre(1, y) + m_mps["PhaseCorrection::C8_3"]->mean() * std::legendre(8, x) * std::legendre(3, y) + m_mps["PhaseCorrection::C9_1"]->mean() * std::legendre(9, x) * std::legendre(1, y);
break;
/*
case 0:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1);
break;
case 1:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1);
break;
case 2:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1);
break;
case 3:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1);
break;
case 4:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1);
break;
case 5:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C0_11"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_9"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_7"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_5"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C4_3"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C5_1"]->mean() * fastLegendre(x,5) * fastLegendre(y,1);
break;
case 6:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C0_11"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C0_13"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_9"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C1_11"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_7"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C2_9"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_5"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C3_7"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C4_3"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C4_5"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C5_1"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C5_3"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C6_1"]->mean() * fastLegendre(x,6) * fastLegendre(y,1);
break;
case 7:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C0_11"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C0_13"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C0_15"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_9"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C1_11"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C1_13"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_7"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C2_9"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C2_11"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_5"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C3_7"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C3_9"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C4_3"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C4_5"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C4_7"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C5_1"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C5_3"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C5_5"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C6_1"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C6_3"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C7_1"]->mean() * fastLegendre(x,7) * fastLegendre(y,1);
break;
case 8:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C0_11"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C0_13"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C0_15"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C0_17"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_9"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C1_11"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C1_13"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C1_15"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_7"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C2_9"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C2_11"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["PhaseCorrection::C2_13"]->mean() * fastLegendre(x,2) * fastLegendre(y,13) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_5"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C3_7"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C3_9"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C3_11"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C4_3"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C4_5"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C4_7"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C4_9"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["PhaseCorrection::C5_1"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C5_3"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C5_5"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C5_7"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["PhaseCorrection::C6_1"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C6_3"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C6_5"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["PhaseCorrection::C7_1"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["PhaseCorrection::C7_3"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["PhaseCorrection::C8_1"]->mean() * fastLegendre(x,8) * fastLegendre(y,1);
break;
case 9:
    return m_mps["PhaseCorrection::C0_1"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C0_3"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C0_5"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C0_7"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C0_9"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C0_11"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C0_13"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C0_15"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C0_17"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["PhaseCorrection::C0_19"]->mean() * fastLegendre(x,0) * fastLegendre(y,19) + m_mps["PhaseCorrection::C1_1"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C1_3"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C1_5"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C1_7"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C1_9"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C1_11"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C1_13"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C1_15"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["PhaseCorrection::C1_17"]->mean() * fastLegendre(x,1) * fastLegendre(y,17) + m_mps["PhaseCorrection::C2_1"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C2_3"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C2_5"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C2_7"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C2_9"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C2_11"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["PhaseCorrection::C2_13"]->mean() * fastLegendre(x,2) * fastLegendre(y,13) + m_mps["PhaseCorrection::C2_15"]->mean() * fastLegendre(x,2) * fastLegendre(y,15) + m_mps["PhaseCorrection::C3_1"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C3_3"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C3_5"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C3_7"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C3_9"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C3_11"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["PhaseCorrection::C3_13"]->mean() * fastLegendre(x,3) * fastLegendre(y,13) + m_mps["PhaseCorrection::C4_1"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C4_3"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C4_5"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C4_7"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C4_9"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["PhaseCorrection::C4_11"]->mean() * fastLegendre(x,4) * fastLegendre(y,11) + m_mps["PhaseCorrection::C5_1"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C5_3"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C5_5"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C5_7"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["PhaseCorrection::C5_9"]->mean() * fastLegendre(x,5) * fastLegendre(y,9) + m_mps["PhaseCorrection::C6_1"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C6_3"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C6_5"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["PhaseCorrection::C6_7"]->mean() * fastLegendre(x,6) * fastLegendre(y,7) + m_mps["PhaseCorrection::C7_1"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["PhaseCorrection::C7_3"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["PhaseCorrection::C7_5"]->mean() * fastLegendre(x,7) * fastLegendre(y,5) + m_mps["PhaseCorrection::C8_1"]->mean() * fastLegendre(x,8) * fastLegendre(y,1) + m_mps["PhaseCorrection::C8_3"]->mean() * fastLegendre(x,8) * fastLegendre(y,3) + m_mps["PhaseCorrection::C9_1"]->mean() * fastLegendre(x,9) * fastLegendre(y,1);
break;
*/

                default:
                    real_t p = 0;
                    for (size_t i=0;i<m_order+1;i++){
                        for (size_t j=0;j<m_order+1-i;j++){
                            //p+=m_mps["PhaseCorrection::C"+std::to_string(i) + "_" + std::to_string(2*j+1)]->mean() * std::legendre(i, x) * std::legendre(2 * j + 1, y);
                            p+=m_mps["PhaseCorrection::C"+std::to_string(i) + "_" + std::to_string(2*j+1)]->mean() * fastLegendre(x,i) * fastLegendre(y,2 * j + 1);
                        }
                    }
                    return p;

                    break;
                }
            }

            


            PhaseCorrection(const MinuitParameterSet& mps) : m_mps(mps),
             m_order(NamedParameter<size_t>("PhaseCorrection::Order", 2)), 
             m_debug(NamedParameter<bool>("PhaseCorrection::debug", true)),
             m_start(NamedParameter<size_t>("PhaseCorrection::Start", 0)),
             m_i(NamedParameter<size_t>("PhaseCorrection::i", 0)),
             m_j(NamedParameter<size_t>("PhaseCorrection::j", 1)),
             m_constrainTransform(NamedParameter<bool>("PhaseCorrection::constrainTransform", false)),
             m_PolyType(NamedParameter<std::string>("PhaseCorrection::PolyType", "antiSym_legendre")) {
                size_t i=0;
            for (size_t j=0;j<m_order+1;j++){
                    for (size_t k=0;k<m_order +1 -j ; k++){
                        auto v= std::vector<size_t>({j, k});
                        auto p = std::pair<size_t, std::vector<size_t> >(i, v);
                        m_ijs.insert(p);
                        if (m_debug) INFO("i = "<<i<<", j = "<<j<<", k = "<<k);
                        i++;
                        
                    }
                }
           }
            PhaseCorrection() = default;
            const std::vector<real_t> getXY(Event event) const {
                auto x = event.s(0,1);
                auto y = event.s(0,2);

                auto z1 = (x + y)/2;
                auto z2 = (y - x)/2;

                real_t m1 = 2.2340742171132946;
                real_t c1 = -3.1171885586526695;

                real_t m2 = 0.8051636393861085;
                real_t c2 = -9.54231895051727e-05;

                auto w1 = m1 * z1 + c1;
                auto w2 = m2 * z2 + c2;


                /*
                 auto x = event.s(0,1);
                auto y = event.s(0,2);
                auto z = event.s(1,2);
                
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
                
                real_t xmin = 0.406004;
                real_t xmax = 2.976556;
                real_t X = (2 * x - xmax - xmin)/(xmax - xmin);
                real_t Y = (2 * y - xmax - xmin)/(xmax - xmin);
                auto zp = (X+Y)/2;
                auto zm = (X-Y)/2;
                */
                return std::vector<real_t>({w1, w2});

            }
            const real_t calcCorr(const Event event)const {
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
                        sum_i += m_mps["PhaseCorrection::C"+std::to_string(i) + std::to_string(j)]->mean()*fastLegendre(zp,i) * fastLegendre(zm,2*j+1);
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


            const Expression fastCorr(const Event event) const {
                auto XY=getXY(event);
                auto x = XY[0];
                auto y = XY[1];
                switch(m_order){
                    case 0:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1);
                    case 1:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1);
                    case 2:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1);
                    case 3:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1);
                    case 4:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1);
                    case 5:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1);
                    case 6:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1);
                    case 7:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C25"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1);
                    case 8:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C08"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C17"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C25"]->mean() * fastLegendre(x,2) * fastLegendre(y,11) + m_mps["PhaseCorrection::C26"]->mean() * fastLegendre(x,2) * fastLegendre(y,13) + m_mps["PhaseCorrection::C30"]->mean() * fastLegendre(x,3) * fastLegendre(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fastLegendre(x,3) * fastLegendre(y,3) + m_mps["PhaseCorrection::C32"]->mean() * fastLegendre(x,3) * fastLegendre(y,5) + m_mps["PhaseCorrection::C33"]->mean() * fastLegendre(x,3) * fastLegendre(y,7) + m_mps["PhaseCorrection::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C35"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C44"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["PhaseCorrection::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C53"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["PhaseCorrection::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C62"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["PhaseCorrection::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["PhaseCorrection::C71"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["PhaseCorrection::C80"]->mean() * fastLegendre(x,8) * fastLegendre(y,1);
                    case 9:
                        return m_mps["PhaseCorrection::C00"]->mean() * fastLegendre(x,0) * fastLegendre(y,1) + m_mps["PhaseCorrection::C01"]->mean() * fastLegendre(x,0) * fastLegendre(y,3) + m_mps["PhaseCorrection::C02"]->mean() * fastLegendre(x,0) * fastLegendre(y,5) + m_mps["PhaseCorrection::C03"]->mean() * fastLegendre(x,0) * fastLegendre(y,7) + m_mps["PhaseCorrection::C04"]->mean() * fastLegendre(x,0) * fastLegendre(y,9) + m_mps["PhaseCorrection::C05"]->mean() * fastLegendre(x,0) * fastLegendre(y,11) + m_mps["PhaseCorrection::C06"]->mean() * fastLegendre(x,0) * fastLegendre(y,13) + m_mps["PhaseCorrection::C07"]->mean() * fastLegendre(x,0) * fastLegendre(y,15) + m_mps["PhaseCorrection::C08"]->mean() * fastLegendre(x,0) * fastLegendre(y,17) + m_mps["PhaseCorrection::C09"]->mean() * fastLegendre(x,0) * fastLegendre(y,19) + m_mps["PhaseCorrection::C10"]->mean() * fastLegendre(x,1) * fastLegendre(y,1) + m_mps["PhaseCorrection::C11"]->mean() * fastLegendre(x,1) * fastLegendre(y,3) + m_mps["PhaseCorrection::C12"]->mean() * fastLegendre(x,1) * fastLegendre(y,5) + m_mps["PhaseCorrection::C13"]->mean() * fastLegendre(x,1) * fastLegendre(y,7) + m_mps["PhaseCorrection::C14"]->mean() * fastLegendre(x,1) * fastLegendre(y,9) + m_mps["PhaseCorrection::C15"]->mean() * fastLegendre(x,1) * fastLegendre(y,11) + m_mps["PhaseCorrection::C16"]->mean() * fastLegendre(x,1) * fastLegendre(y,13) + m_mps["PhaseCorrection::C17"]->mean() * fastLegendre(x,1) * fastLegendre(y,15) + m_mps["PhaseCorrection::C18"]->mean() * fastLegendre(x,1) * fastLegendre(y,17) + m_mps["PhaseCorrection::C20"]->mean() * fastLegendre(x,2) * fastLegendre(y,1) + m_mps["PhaseCorrection::C21"]->mean() * fastLegendre(x,2) * fastLegendre(y,3) + m_mps["PhaseCorrection::C22"]->mean() * fastLegendre(x,2) * fastLegendre(y,5) + m_mps["PhaseCorrection::C23"]->mean() * fastLegendre(x,2) * fastLegendre(y,7) + m_mps["PhaseCorrection::C24"]->mean() * fastLegendre(x,2) * fastLegendre(y,9) + m_mps["PhaseCorrection::C25"]->mean() * fcn::pow(x,2) * fcn::pow(y,11) + m_mps["PhaseCorrection::C26"]->mean() * fcn::pow(x,2) * fcn::pow(y,13) + m_mps["PhaseCorrection::C27"]->mean() * fcn::pow(x,2) * fcn::pow(y,15) + m_mps["PhaseCorrection::C30"]->mean() * fcn::pow(x,3) * fcn::pow(y,1) + m_mps["PhaseCorrection::C31"]->mean() * fcn::pow(x,3) * fcn::pow(y,3) + m_mps["PhaseCorrection::C32"]->mean() * fcn::pow(x,3) * fcn::pow(y,5) + m_mps["PhaseCorrection::C33"]->mean() * fcn::pow(x,3) * fcn::pow(y,7) + m_mps["PhaseCorrection::C34"]->mean() * fastLegendre(x,3) * fastLegendre(y,9) + m_mps["PhaseCorrection::C35"]->mean() * fastLegendre(x,3) * fastLegendre(y,11) + m_mps["PhaseCorrection::C36"]->mean() * fastLegendre(x,3) * fastLegendre(y,13) + m_mps["PhaseCorrection::C40"]->mean() * fastLegendre(x,4) * fastLegendre(y,1) + m_mps["PhaseCorrection::C41"]->mean() * fastLegendre(x,4) * fastLegendre(y,3) + m_mps["PhaseCorrection::C42"]->mean() * fastLegendre(x,4) * fastLegendre(y,5) + m_mps["PhaseCorrection::C43"]->mean() * fastLegendre(x,4) * fastLegendre(y,7) + m_mps["PhaseCorrection::C44"]->mean() * fastLegendre(x,4) * fastLegendre(y,9) + m_mps["PhaseCorrection::C45"]->mean() * fastLegendre(x,4) * fastLegendre(y,11) + m_mps["PhaseCorrection::C50"]->mean() * fastLegendre(x,5) * fastLegendre(y,1) + m_mps["PhaseCorrection::C51"]->mean() * fastLegendre(x,5) * fastLegendre(y,3) + m_mps["PhaseCorrection::C52"]->mean() * fastLegendre(x,5) * fastLegendre(y,5) + m_mps["PhaseCorrection::C53"]->mean() * fastLegendre(x,5) * fastLegendre(y,7) + m_mps["PhaseCorrection::C54"]->mean() * fastLegendre(x,5) * fastLegendre(y,9) + m_mps["PhaseCorrection::C60"]->mean() * fastLegendre(x,6) * fastLegendre(y,1) + m_mps["PhaseCorrection::C61"]->mean() * fastLegendre(x,6) * fastLegendre(y,3) + m_mps["PhaseCorrection::C62"]->mean() * fastLegendre(x,6) * fastLegendre(y,5) + m_mps["PhaseCorrection::C63"]->mean() * fastLegendre(x,6) * fastLegendre(y,7) + m_mps["PhaseCorrection::C70"]->mean() * fastLegendre(x,7) * fastLegendre(y,1) + m_mps["PhaseCorrection::C71"]->mean() * fastLegendre(x,7) * fastLegendre(y,3) + m_mps["PhaseCorrection::C72"]->mean() * fastLegendre(x,7) * fastLegendre(y,5) + m_mps["PhaseCorrection::C80"]->mean() * fastLegendre(x,8) * fastLegendre(y,1) + m_mps["PhaseCorrection::C81"]->mean() * fastLegendre(x,8) * fastLegendre(y,3) + m_mps["PhaseCorrection::C90"]->mean() * fastLegendre(x,9) * fastLegendre(y,1);
             
                    default:
                        return calcCorr(event);
                }



            }


            const real_t calcCorrXY(double x, double y) const {
                double p = 0;
                  auto z1 = (x + y)/2;
                    auto z2 = (y - x)/2;

    real_t m1 = 2.2340742171132946;
                    real_t c1 = -3.1171885586526695;

                    real_t m2 = 0.8051636393861085;
                    real_t c2 = -9.54231895051727e-05;

                    auto w1 = m1 * z1 + c1;
                    auto w2 = m2 * z2 + c2;


                if (m_PolyType=="Gaussian"){
                    real_t sc = m_mps["PhaseCorrection::GaussScale"]->mean();
                    real_t muX = m_mps["PhaseCorrection::GaussMuX"]->mean();
                    real_t muY = m_mps["PhaseCorrection::GaussMuY"]->mean();
                    real_t sigmaX = m_mps["PhaseCorrection::GaussSigmaX"]->mean();
                    real_t sigmaY = m_mps["PhaseCorrection::GaussSigmaY"]->mean();
                    //real_t linWm = m_mps["PhaseCorrection::GaussLinearW-"]->mean();
                    //real_t quadWm = m_mps["PhaseCorrection::GaussQuadraticW-"]->mean();
                    real_t erfFactor = m_mps["PhaseCorrection::GaussErfFactor"]->mean();
                    int sign = z2/abs(z2);
                    real_t mod_erf = std::erf(w2/erfFactor);
                    if (sign == 1){
                    //real_t gauss = (linWm * wm + sign * quadWm * pow(w2, 2)) * sc/(2 * M_PI * sigmaX * sigmaY) * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY));
                    //real_t gauss = sign * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;
                    //real_t erf = 1 - pow(1 + 0.278393 * abs(w2) + 0.230389 * pow(abs(w2), 2) + 0.000972 * pow(abs(w2), 3) + 0.078108 * pow(abs(w2), 4), -4);
                    
                    //real_t gauss =  (linWm * w2 + sign * quadWm * pow(w2, 2)) * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;
                    real_t gauss =  mod_erf * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;


                    return gauss;

                    }
                    else{
                    //real_t gauss = sign * sc/(2 * M_PI * sigmaX * sigmaY) * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    //real_t gauss = sign * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    //real_t gauss = (linWm * w2 + sign * quadWm * pow(w2, 2)) * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    real_t gauss = mod_erf * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));


                    return gauss;
 
                    }

                }
/*
                for (size_t i=m_start;i<m_order+1 ; i++){
                  
                 
                    p+=calcPoly(w1,w2,i);
                    if (m_debug) INFO("f"<<i<<" = "<<p);

                }
            */

                p = Legendre2D(w1, w2);
                if (m_constrainTransform){
                    if (p>M_PI) p = p - 2 * M_PI;
                    if (p<-M_PI) p = p + 2 * M_PI;
                }
                return p;
                

            }

            const real_t calcCorrL(const Event event) const {
                //auto XY=getXY(event);
                //auto x = XY[0]* 10;
                //auto y = XY[1];
                auto x = event.s(0,1);
                auto y = event.s(0,2);

                auto z1 = (x + y)/2;
                auto z2 = (y - x)/2;

                real_t m1 = 2.2340742171132946;
                real_t c1 = -3.1171885586526695;

                real_t m2 = 0.8051636393861085;
                real_t c2 = -9.54231895051727e-05;

                auto w1 = m1 * z1 + c1;
                auto w2 = m2 * z2 + c2;

                real_t p = 0;
                
               
                /*
                int i=m_i;
                int j=m_j;

                

                return m_mps["PhaseCorrection::C" + std::to_string(i) + "_"  + std::to_string(j)]->mean() * fastLegendre(w1, i) * fastLegendre(w2, j);
                */
                
                if (m_PolyType=="Gaussian"){
                    real_t sc = m_mps["PhaseCorrection::GaussScale"]->mean();
                    real_t muX = m_mps["PhaseCorrection::GaussMuX"]->mean();
                    real_t muY = m_mps["PhaseCorrection::GaussMuY"]->mean();
                    real_t sigmaX = m_mps["PhaseCorrection::GaussSigmaX"]->mean();
                    real_t sigmaY = m_mps["PhaseCorrection::GaussSigmaY"]->mean();
                    //real_t linWm = m_mps["PhaseCorrection::GaussLinearW-"]->mean();
                    //real_t quadWm = m_mps["PhaseCorrection::GaussQuadraticW-"]->mean();
                    real_t erfFactor = m_mps["PhaseCorrection::GaussErfFactor"]->mean();
                    int sign = z2/abs(z2);
                    real_t mod_erf = std::erf(w2/erfFactor);
                    if (sign == 1){
                    //real_t gauss = (linWm * wm + sign * quadWm * pow(w2, 2)) * sc/(2 * M_PI * sigmaX * sigmaY) * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY));
                    //real_t gauss = sign * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;
                    //real_t erf = 1 - pow(1 + 0.278393 * abs(w2) + 0.230389 * pow(abs(w2), 2) + 0.000972 * pow(abs(w2), 3) + 0.078108 * pow(abs(w2), 4), -4);
                    
                    //real_t gauss =  (linWm * w2 + sign * quadWm * pow(w2, 2)) * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;
                    real_t gauss =  mod_erf * sc * exp ( -pow( x - muX, 2)/(2*sigmaX) - pow(y - muY, 2)/(2 * sigmaY)) ;


                    return gauss;

                    }
                    else{
                    //real_t gauss = sign * sc/(2 * M_PI * sigmaX * sigmaY) * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    //real_t gauss = sign * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    //real_t gauss = (linWm * w2 + sign * quadWm * pow(w2, 2)) * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));
                    real_t gauss = mod_erf * sc * exp ( -pow( y - muX, 2)/(2*sigmaX) - pow(x - muY, 2)/(2 * sigmaY));


                    return gauss;
 
                    }

                }
                /*
                for (size_t i=m_start;i<m_order+1 ; i++){
                   
                    p+=calcPoly(w1,w2,i);
                    if (m_debug) INFO("f"<<i<<" = "<<p);

                }
                */
                p = Legendre2D(w1, w2);
                return p;
                
                
                //return m_mps["PhaseCorrection::C" + std::to_string(i) + "_"  + std::to_string(j)]->mean() * std::pow(x, i) * std::pow(y, j);

            }

            const real_t est_errCorrL(const Event event, TMatrixTSym<double>* pcov,  bool useCov=false , size_t startFrom=0, size_t until=0) const {
                real_t result = 0;
                
                auto x = event.s(0,1);
                auto y = event.s(0,2);

                auto z1 = (x + y)/2;
                auto z2 = (y - x)/2;

                real_t m1 = 2.2340742171132946;
                real_t c1 = -3.1171885586526695;

                real_t m2 = 0.8051636393861085;
                real_t c2 = -9.54231895051727e-05;

                auto w1 = m1 * z1 + c1;
                auto w2 = m2 * z2 + c2;

                if (m_PolyType=="Gaussian"){
                    return 0;
                }



                
              //  INFO("At start of err_estCorrL"); 
 

               auto ijs = m_ijs;
                        

                if (useCov){

                    auto cov = (*pcov);

                for (size_t i=startFrom; i<cov.GetNrows() - until;i++){
                    for (size_t j=startFrom; j<cov.GetNcols() - until;j++){
                        size_t i1 = ijs[i][0];
                        size_t j1 = ijs[i][1];
                        size_t i2 = ijs[j][0];
                        size_t j2 = ijs[j][1];
                        //result += cov[i][j] * Poly1D(x, i) * Poly1D(y, 2*j+1);
                        result += cov[i][j] * Poly1D(w1, i1) * Poly1D(w2, 2*j1+1) * Poly1D(w1, i2) * Poly1D(w2, 2*j2+1);
                    }
                }


                }
                else{
                    /*
                    int i=m_i;
                    int j=m_j;
                    result += std::norm(m_mps["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(j)]->err()) * Poly1D(x, i) * Poly1D(y, j) * Poly1D(x, i) * Poly1D(y, j) ;
                    */

                    
                    for (size_t i=0; i < m_order; i++){
                        for (size_t j=0;j< m_order -i;j++){
                            result += std::norm( m_mps["PhaseCorrection::C" + std::to_string(i) + "_" + std::to_string(2 * j +1)]->err() * Poly1D(w1, i) * Poly1D(w2, 2*j+1) );
                        }
                    }
                    
                    

                }
                return sqrt(result);
            }

            const real_t est_errCov(const Event event1, const Event event2, TMatrixTSym<double> cov = TMatrixTSym<double>(), size_t startFrom=0, size_t until=0) const {
                real_t result = 0;
                auto XY1 = getXY(event1);
                auto x1 = XY1[0];
                auto y1 = XY1[1];

                auto XY2 = getXY(event2);
                auto x2 = XY2[0];
                auto y2 = XY2[1];



                
                INFO("At start of err_estCorrL"); 
 

                int nE = cov.GetNoElements();
                INFO("nE = "<<nE);
              


                auto ijs = m_ijs;
                        

                if (nE>0){

                for (size_t i=startFrom; i<cov.GetNrows() - until;i++){
                    for (size_t j=startFrom; j<cov.GetNcols() - until;j++){
             
             
             
                        size_t i1 = ijs[i][0];
             
                        size_t j1 = ijs[i][1];
                        size_t i2 = ijs[j][0];
                        size_t j2 = ijs[j][1];


                     


                        //result += cov[i][j] * Poly1D(x, i) * Poly1D(y, 2*j+1);
                        result += cov[i][j] * Poly1D(x1, i1) * Poly1D(y1, 2*j1+1) * Poly1D(x2, i2) * Poly1D(y2, 2*j2+1);

                    }
                }
                }
                else{
                    int i=m_i;
                    int j=m_j;
                    result += std::norm(m_mps["PhaseCorrecion::C" + std::to_string(i) + "_" + std::to_string(j)]->err()) * Poly1D(x1, i) * Poly1D(y1, j) * Poly1D(x2, i) * Poly1D(y2, j) ;
                    /*
                    for (size_t i=0; i < m_order; i++){
                        for (size_t j=0; m_order -i;j++){
                            result += std::norm(m_mps["PhaseCorrecion::C" + std::to_string(i) + "_" + std::to_string(2 * j +1)]->err()) * Poly1D(x1, i) * Poly1D(y1, 2*j+1) * Poly1D(x2, i) * Poly1D(y2, 2*j+1) ;

                        }
                    }
                    */

                }
                return result;
            }

            std::function<real_t(real_t,real_t,size_t,MinuitParameterSet)> makefunc(){
                size_t order = m_order;
                std::function<real_t(real_t,real_t,int,MinuitParameterSet)> f = 
                    [](real_t x, real_t y,size_t order, MinuitParameterSet mps){
                        real_t corr=0; 
                        
                        for (auto i=0; i < order+1; i++){
                        real_t sum_i=0;
                            for (auto j=0; j<order+1-i; j++){
                            real_t Cij = mps["PhaseCorrection::C"+std::to_string(i) + std::to_string(j)]->mean();
                            
                            //sum_i += Cij*fastLegendre(x,i) * fastLegendre(y,2*j+1);
                            sum_i += Cij*std::pow(x,i)* std::pow(y,2*j+1);
                            }
                        corr += sum_i;
                        }
                        
                    return corr;
                };
                return f;
            }
            

        void setEvents(EventList& list){
            m_events = &list;
        }
        
        void prepareCache(){
            INFO("Prepare cache for order "<<m_order);

            int id=0;
            for (auto& evt : *m_events){
                auto XY=getXY(evt);
                auto x = XY[0];
                auto y = XY[1];
                if (id%10000==0)std::cout<<"\rAt "<<id<<" event out of"<<m_events->size()<<std::flush;
                id++;
                std::map<std::string, real_t> map;
                for (auto i=0;i<m_order+1;i++){
                    for (auto j=0;j<m_order+1-i;j++){
                        std::stringstream name_ss;
                        name_ss<<"PhaseCorrection::C"<<i<<"_"<<2*j+1;
                        std::string name_ij = name_ss.str();
                       
                        real_t value_ij = fastLegendre(x, i) * fastLegendre(y, 2*j+1);
                       
                        map.insert(std::pair<std::string, real_t> (name_ij, value_ij));
                    }
                }
              
                

                //m_cache.insert(std::pair<const real_t*, const std::map<std::string, real_t> > (evt.address(), map));
                m_cache.push_back(map);
            }
            INFO("Prepared Cache");

        }


        //const std::map<const real_t*, const std::map<std::string, real_t> > getCache(){
        const std::vector<std::map<std::string, real_t> > getCache(){
            return m_cache;
        };

        const real_t getValCache(int i)const{
            real_t v=0;
            auto m = m_cache[i];
            for (auto& [key, value] : m){
                v+=m_mps[key]->mean() * m.at(key);
            }
            return v;

        }
        void setOrder(size_t order){
            m_order = order;
        }
        

            

            
        private:
            MinuitParameterSet m_mps;
            size_t m_order;
            size_t m_start;
            size_t m_i;
            size_t m_j;
            bool m_debug;
            bool m_constrainTransform;
	    std::string m_PolyType;
            EventList* m_events = {nullptr};
            //std::map< const real_t*, const std::map<std::string, real_t> > m_cache = {};
            std::vector< std::map<std::string, real_t> > m_cache = {};
            std::map<size_t , std::vector<size_t> > m_ijs;

          

    };
}
#endif
