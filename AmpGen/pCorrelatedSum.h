#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"

#ifndef PCORRELATEDSUM
#define PCORRELATEDSUM

/* 
   pCorrelatedSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
   This object is designed to take the correlated decays:
   P --> D1 D2
   D1 --> X1 and D2 --> X2 or D2 --> X1 and D1 --> X2
   which will give an amplitude in the X1 phasespace (x1) and the X2 phasespace (x2)

   A(x1, x2) = A(x1) B(x2) - C(x1) D(x2) [*]

   where A = Amp(D1-->X1), B = Amp(D2 --> X2), C = Amp(D2 --> X1) and D = Amp(D1 --> X2)
   the corresponding probability density is thus

   P(x1,x2) = |A(x1)|^2 |B(x2)|^2 + |C(x1)|^2 |D(x2)|^2 - 2 Re(A(x1) B*(x2) C(x1) D*(x2))
   Therefore to normalise this distribution we need to integrate P(x1, x2):

   N = int_x1 int_x2 P(x1, x2) dx1 dx2 
   = NA NB + NC ND - 2 Re(NAC NBD)

   where N{FG} = int_x F(x) G*(x) dx
   and N{FF} = N{F} = int_x |F(x)|^2 dx

   These "N" factors are calculated by summing individual components of the amplitude so 
   N{F} = sum_i sum_j g_i g*_j NA_ij
   where the indicies are for each amplitude A_i(x) and g_i is the coupling to that component



*/


namespace AmpGen { 

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


class pCorrelatedSum {
    public:
      //Takes amplitudes - 
        pCorrelatedSum();
        pCorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps);
        virtual ~pCorrelatedSum()=default;
        real_t operator()( const Event& event1, const Event& event2) const { return prob(event1, event2); }
        real_t prob(const Event& event1, const Event& event2) const {
        double P = std::norm(getVal(event1, event2))/m_norm;  
        double A2 = std::norm(m_A.getVal(event1))/m_norm;
        double B2 = std::norm(m_B.getVal(event2));
        double C2 = std::norm(m_C.getVal(event1))/m_norm;
        double D2 = std::norm(m_D.getVal(event2));
        complex_t AC = m_A.getVal(event1) * std::conj(m_C.getVal(event1))/m_norm;
        complex_t BD = m_B.getVal(event2) * std::conj(m_D.getVal(event2));
        auto i = Constant(0,1);
        complex_t eif = exp(i() * correction(event1));
        if (m_sameTag){
            eif *= exp(-i() * correction(event2));
        }
        double inter = -2*std::real(AC * BD);
        double corrected_inter = -2*std::real(AC * BD * eif);
        if (m_debug){
            INFO("|AB-CD|^2 = "<<P);
            INFO("A^2 B^2 + C^2 D^2 - 2Re(AC*BD) = "<<A2 *B2+C2*D2+inter);
            INFO("A^3 B^2 + C^2 D^2 - 2Re(AC*BD*eif) = "<<A2 *B2+C2*D2+corrected_inter);
            INFO("|A|^2 = "<<A2);
            INFO("|B|^2 = "<<B2);
            INFO("|C|^2 = "<<C2);
            INFO("|D|^2 = "<<D2);
            INFO("-2ReAC*BD*eif = "<<inter);
            INFO("Strong phase (AC) = "<<std::imag(std::log(AC/std::abs(AC))));
            INFO("Correction = "<<std::imag(std::log(eif)));
            INFO("Corrected Strong phase (AC) = "<<std::imag(std::log(AC*eif/std::abs(AC))));
            INFO("Norm = "<<m_norm);

        }
        return P;
        }

    double getC(int i, int j)const {
        std::string key = "pCorrelatedSum::C"+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = m_mps[key]->mean(); 
        return val;
    }


    real_t prob_unnormalised(const Event& event1, const Event& event2) const {return std::norm(getVal(event1, event2));}
    void prepare();
    void reset(bool resetEvents);
    double slowNorm();
    void setEvents(EventList& list1, EventList& list2);
    void setMC(EventList& list1, EventList& list2);

    void updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iB,
        const std::vector<size_t>& iC, const std::vector<size_t>& iD);

    void debugNorm();
    void debug(const Event& event1, const Event& event2) const; 
    complex_t getVal(const Event& event1, const Event& event2) const;
    complex_t getValNoCache(const Event& event1, const Event& event2) const;
    complex_t getValNoCache(const Event& event1, const Event& event2, const size_t& offset) const;
    real_t size()const { return m_A.size() + m_B.size() + m_C.size() + m_D.size();}
    std::vector<std::vector<FitFraction> > fitFractions(const LinearErrorPropagator& linProp);
    real_t norm()  const;
    double probA(const Event& event){
        double prob = m_A.prob(event);
        return prob;
    }
    std::vector<complex_t> getVals(const Event& event1, const Event& event2) const {
        complex_t A = m_A.getVal(event1);
        complex_t B = m_B.getVal(event2);
        complex_t C = m_C.getVal(event1);
        complex_t D = m_D.getVal(event2);
        complex_t ABCD = getVal(event1, event2);
        complex_t corr = correction(event1);
        std::vector<complex_t> vals = {A,B,C,D,ABCD, corr};
        return vals;
    }
    complex_t correction(const Event& event) const {
        Expression corr = 0;
        auto x = event.s(0,1);
        auto y = event.s(0,2);
        auto z = event.s(1,2);
        auto mp = sqrt(event.s(1,1));
        auto mm = sqrt(event.s(2,2));
        auto mK = sqrt(event.s(0,0));
        auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
        Expression xmin = pow(mp + mK, 2);
        Expression xmax = pow(mD - mm, 2);
        Expression x0 = (xmax + xmin)/2;
        Expression ymin = pow(mp + mK, 2);
        Expression ymax = pow(mD - mp, 2);
        Expression y0 = (ymax + ymin)/2;
        Expression X = x/x0 - 1;
        Expression Y = y/y0 - 1;
        //INFO("Order = "<<m_order);
        for (int i=0; i < m_order+1; i++){
            Expression sum_i=0;
            for (unsigned int j=0; j<m_order+1-i; j++){
                double Cij = getC(i,j);
                if (m_pdebug){
                    INFO("x = "<<x);
                    INFO("y = "<<y);
                    INFO("Type = "<<m_polyType) ;
                    INFO("C"<<i<<j<<" = "<<Cij) ;
                }
                if (m_polyType=="simple"){
                    sum_i = sum_i +  Cij * pow(X(), i) * pow(Y(), j);
                }
                else if (m_polyType=="chebychev"){
                sum_i = sum_i + Cij * chebychev(X(), i) * chebychev(Y(), j); 
                }
                else if (m_polyType=="legendre"){
                    sum_i = sum_i + Cij * legendre(X(), i) * legendre(Y(), j);
                }
                else if (m_polyType=="laguerre"){
                    sum_i = sum_i + Cij * laguerre(X(), i) * laguerre(Y(), j);
                }
                if (m_pdebug){
                    INFO("sum_"<<i<<" = "<<sum_i());
                }

            }
            corr += sum_i;
        }
        //corr = Constant(0,1) * corr;
        //complex_t val = exp(corr());
        if (m_pdebug) INFO("correction = "<<corr());
        return corr();
      }
      //real_t norm(const Bilinears& norms) const; 
      double m_inter = 0;

    double getCParam(int i, int j) const{
        std::string key = "MyPoly::C"+std::to_string(i)+std::to_string(j);
        bool present=false;
        for (auto param : m_mps){
            if (param->name()==key){
                present=true;
                return param->mean();
            }
        }
        return 0;
    }
   
    protected:
        double  m_norm  =    {0};
        Bilinears m_normalisationsAC;
        Bilinears m_normalisationsBD;
        Bilinears m_normalisationsAA;
        Bilinears m_normalisationsBB;
        Bilinears m_normalisationsCC;
        Bilinears m_normalisationsDD;
        MinuitParameterSet m_mps;
        CoherentSum  m_A;
        CoherentSum  m_B;
        CoherentSum  m_C;
        CoherentSum  m_D;
        EventList* m_events1 = nullptr;
        EventList* m_events2 = nullptr;
        Integrator<10> m_integratorAC;
        Integrator<10> m_integratorBD;
        std::vector<std::vector<MinuitParameter*> > m_Cparams;
        std::vector<std::vector<double> > m_poly;
        std::vector<FitFraction> outputFractionsA;
        std::vector<FitFraction> outputFractionsB; 
        std::vector<FitFraction> outputFractionsC;
        std::vector<FitFraction> outputFractionsD;
        Integrator<10> m_integratorAA;
        Integrator<10> m_integratorBB;
        Integrator<10> m_integratorCC;
        Integrator<10> m_integratorDD;
        bool m_coherentIntegral;
        bool m_coherentIntegralA;
        bool m_coherentIntegralB;
        bool m_coherentIntegralC;
        bool m_coherentIntegralD;
        bool m_debug;
        int m_debugFreq;
        bool m_pdebug;
        bool m_pNorm;
        bool m_updateNorms;
        int m_order;
        std::string m_polyType;
        size_t m_prepareCalls = 0;
        double m_Anorm;
        double m_Bnorm;
        double m_Cnorm;
        double m_Dnorm;
        bool m_sameTag;
        bool m_flat;
  };
}
#endif
