#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"

#ifndef PCOHERENTSUM
#define PCOHERENTSUM

/* 
   pCoherentSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
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

Expression CPSinPoly2(Expression x, Expression y, unsigned int mu, unsigned int lambda, int CP){
  auto M00 = fcn::sin(M_PI * lambda * x);
  auto M01 = fcn::sin(M_PI * lambda * y);
  auto M10 = fcn::sin(M_PI * mu * x);
  auto M11 = fcn::sin(M_PI * mu * y);
  return M00 * M11 + CP*M10*M01;
}

Expression chebychev2(Expression x, unsigned int order){
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
        output = 2 * x * chebychev2(x, order -1) - chebychev2(x, order - 2);
    }
    return output;
}


Expression legendre2(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x;
    }
    else {
        output = (2 * order - 1)/order * x * legendre2(x, order - 1) - (order - 1)/order * legendre2(x, order - 2);
    }
    return output;
}
    
Expression bessel2(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x + 1;
    }
    else{
        output = (2 * order - 1) * x * bessel2(x, order - 1) + bessel2(x, order - 2);
        
    }
    return output;
}

Expression laguerre2(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output= 1 - x;
    }
    else{
        output = ( (2 * order  - 1 - x) * laguerre2(x, order - 1) - (order - 1) * laguerre2(x, order - 2))/order;
    }
    return output;
}


class pCoherentSum {
    public:
      //Takes amplitudes - 
        pCoherentSum();
        pCoherentSum(const EventType& type1, const MinuitParameterSet& mps, std::string SFType="Psi3770", int gamSign=1, bool useXY=false);
        virtual ~pCoherentSum()=default;


        real_t operator()( const Event& event1) const { return prob(event1); }
        real_t prob(const Event& event1) const {
        double P = std::norm(getVal(event1))/m_norm;  
        double A2 = std::norm(m_A.getVal(event1))/m_norm;

        double C2 = std::norm(m_C.getVal(event1))/m_norm;

        complex_t AC = m_A.getVal(event1) * std::conj(m_C.getVal(event1))/m_norm;

        auto i = Constant(0,1);
        complex_t eif = exp(i() * correction(event1));
        
        double inter = -2*std::real(AC);
        double corrected_inter = -2*std::real(AC * eif);
        if (m_debug){
            INFO("|AB-CD|^2 = "<<P);
            INFO("A^2  + C^2  - 2Re(AC) = "<<A2+C2+inter);
            INFO("A^2 + C^2 - 2Re(AC*eif) = "<<A2 +C2+corrected_inter);
            INFO("|A|^2 = "<<A2);

            INFO("|C|^2 = "<<C2);

            INFO("-2ReAC*eif = "<<inter);
            INFO("Strong phase (AC) = "<<std::imag(std::log(AC/std::abs(AC))));
            INFO("Correction = "<<std::imag(std::log(eif)));
            INFO("Corrected Strong phase (AC) = "<<std::imag(std::log(AC*eif/std::abs(AC))));
            INFO("Norm = "<<m_norm);

        }
        return P;
        }

    double getC(int i, int j, std::string pref="")const {
        std::string key = "pCorrelatedSum::C"+pref+std::to_string(i)+std::to_string(j);
        double val=0;
        //if (m_debug) INFO(m_mps[key]->name()<<" = "<<m_mps[key]->mean());
        val = m_mps[key]->mean(); 
        return val;
    }

    complex_t getSumFactor()const {
        if (m_SFType=="Psi3770"){
            return -1;
        }
        else{

            if (m_useXY){
        if (m_gamSign==1){
        double xp = m_mps["pCoherentSum::x+"]->mean();
        double yp = m_mps["pCoherentSum::y+"]->mean();
       

    
 

        auto sumFactor =  xp + Constant(0, 1)() * yp;
//        auto sumFactor = xp + std::complex<double>(0, 1)* yp ;//+ Constant(0, 1) * yp;


            return sumFactor;
        }
        else {
        double xm = m_mps["pCoherentSum::x-"]->mean();
        double ym = m_mps["pCoherentSum::y-"]->mean();
 
        complex_t sumFactor = xm + Constant(0, 1)() * ym;
            return sumFactor;
        }





            }
else{

        double r = 0;
        double d = 0;
        double g = 0;
        std::stringstream ss_key_r;
        std::stringstream ss_key_d;
        std::stringstream ss_key_g;
        ss_key_r<<"pCoherentSum::"<<m_SFType<<"r";
        ss_key_d<<"pCoherentSum::"<<m_SFType<<"d";
        ss_key_g<<"pCoherentSum::gamma";

        std::string key_r = ss_key_r.str();
        std::string key_d = ss_key_d.str();
        std::string key_g = ss_key_g.str();

        if (m_debug) INFO("Key for rB = "<<key_r);
        if (m_debug) INFO("Key for deltaB = "<<key_d);
        if (m_debug) INFO("Key for gamma = "<<key_g);

        auto mps_r = m_mps.find(key_r);
        auto mps_d = m_mps.find(key_d);
        auto mps_g = m_mps.find(key_g);




            r = mps_r->mean();
        if (m_debug) INFO("rB = "<<r);




            d = mps_d->mean();

        if (m_debug) INFO("deltaB = "<<d);



            g = mps_g->mean() * m_gamSign;

        if (m_debug) INFO("gamma = "<<g);

        complex_t sumFactor = exp(Constant(0, 1)() * (g + d)) * r;
        return sumFactor;
        }
        }
        
    }
    
    

    real_t prob_unnormalised(const Event& event1) const {return std::norm(getVal(event1));}
    void prepare();
    void reset(bool resetEvents);
    double slowNorm();
    real_t testnorm() const;
    void setEvents(EventList& list1);
    void setMC(EventList& list1);

    void updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iC);

    void debugNorm();
    void debug(const Event& event1) const; 
    complex_t getVal(const Event& event1) const;
    complex_t getValNoCache(const Event& event1) const;
    complex_t getValNoCache(const Event& event1, const size_t& offset) const;
    real_t size()const { return m_A.size() + m_C.size() ;}
    std::vector<std::vector<FitFraction> > fitFractions(const LinearErrorPropagator& linProp);
    real_t norm()  const;
    double probA(const Event& event){
        double prob = m_A.prob(event);
        return prob;
    }
    std::vector<complex_t> getVals(const Event& event1) const {
        complex_t A = m_A.getVal(event1);
    
        complex_t C = m_C.getVal(event1);
    
        complex_t AC = getVal(event1);
        complex_t corr = correction(event1);
        std::vector<complex_t> vals = {A,C,AC, corr};
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
        Expression X = (2 * x - xmax - xmin)/(xmax - xmin);
        Expression Y = (2 * y - ymax - ymin)/(ymax - ymin);
        
        //INFO("Order = "<<m_order);
        for (auto i=0; i < m_order+1; i++){
            Expression sum_i=0;
            for (auto j=0; j<m_order+1-i; j++){
                if (m_polyType=="CPSinPoly"){
                  double Cpij = getC(i,j,"P");
                  double Cmij = getC(i,j,"M");
                  sum_i = sum_i + Cpij * CPSinPoly2(X(), Y(), i, j, 1) + Cmij * CPSinPoly2(X(), Y(), i, j, -1);
                }
                else if (m_polyType=="CP_legendre"){
                   double Cpij = getC(i,j,"P");
                   double Cmij = getC(i,j,"M");
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i + Cpij* legendre2(zp, i) * legendre2(zm, 2*j) + Cmij * legendre2(zp, i) * legendre2(zm, 2*j+1);
                   

                }
                else if (m_polyType=="antiSym_legendre"){
                   double Cij = getC(i,j);

                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * legendre2(zp, i) * legendre2(zm, 2*j+1);
                   

                }
                else if (m_polyType=="antiSym_simple"){
                   double Cij = getC(i,j);
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * std::pow(zp, i) * std::pow(zm, 2*j+1);
                   

                }
                else if (m_polyType=="antiSym_chebyshev"){
                   double Cij = getC(i,j);
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * chebychev2(zp, i) * chebychev2(zm, 2*j+1);
                   

                }
                else if (m_polyType=="Sym_legendre"){
                   double Cij = getC(i,j);
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * legendre2(zp, i) * legendre2(zm, 2*j);
                   

                }



                else{
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
                sum_i = sum_i + Cij * chebychev2(X(), i) * chebychev2(Y(), j); 
                }
                else if (m_polyType=="legendre"){
                    sum_i = sum_i + Cij * legendre2(X(), i) * legendre2(Y(), j);
                }
                else if (m_polyType=="laguerre"){
                    sum_i = sum_i + Cij * laguerre2(X(), i) * laguerre2(Y(), j);
                }
                else if (m_polyType=="bessel"){
                    sum_i = sum_i + Cij * bessel2(X(), i) * bessel2(Y(), j);
                }
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

      complex_t A(Event& evt){
          return m_A.getVal(evt);
      }
      complex_t C(Event& evt){
          return m_C.getVal(evt);
      }

      real_t getFastNorm(){

        if (m_debug) INFO("Getting the value for the normalised pdf");
        if (m_debug) INFO("Get Matrix Elements for A,B,C,D");
        complex_t sumFactor = getSumFactor();  
        auto sum_amps = []( const Bilinears& bl, const auto& mA, const auto& mB )
        {
            complex_t v; 
            for (size_t i=0; i< mA.size(); i++){
            for (size_t j=0; j< mB.size(); j++){
                v += bl.get(i, j) * mA[i].coefficient * std::conj(mB[j].coefficient);;
            }
            }
            return v; 
        };
        /*

            */
        real_t nA = m_A.norm();

        real_t nC = m_C.norm();


        complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );
        complex_t mix = nAC * std::conj(sumFactor);
        real_t intTerm = 2 * mix.real();


        double xp = m_mps["pCoherentSum::x+"]->mean();
        double yp = m_mps["pCoherentSum::y+"]->mean();

        
        real_t N = nA + nC * std::abs(sumFactor) + intTerm;

        
        

        return N;


      }
   
    protected:
        double  m_norm  =    {0};
        Bilinears m_normalisationsAC;

        Bilinears m_normalisationsAA;

        Bilinears m_normalisationsCC;

        MinuitParameterSet m_mps;
        std::string m_SFType;
        CoherentSum  m_A;

        CoherentSum  m_C;

        EventList* m_events1 = nullptr;

        Integrator<10> m_integratorAC;

        std::vector<std::vector<MinuitParameter*> > m_Cparams;
        std::vector<std::vector<double> > m_poly;
        std::vector<FitFraction> outputFractionsA;

        std::vector<FitFraction> outputFractionsC;

        Integrator<10> m_integratorAA;

        Integrator<10> m_integratorCC;

        bool m_coherentIntegral;
        bool m_coherentIntegralA;

        bool m_coherentIntegralC;

        bool m_debug;
        int m_debugFreq;
        bool m_pdebug;
        bool m_pNorm;
        bool m_fastNorm;
        bool m_updateNorms;
        int m_order;
        std::string m_polyType;
        size_t m_prepareCalls = 0;
        double m_Anorm;

        double m_Cnorm;

        int m_gamSign;

        bool m_sameTag;
        bool m_flat;
        bool m_useXY;
  };
}
#endif

