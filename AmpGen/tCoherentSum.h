#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MinuitParameterSet.h"

#ifndef TCOHERENTSUM
#define TCOHERENTSUM

/* 
   tCoherentSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
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

Expression CPSinPolyt(Expression x, Expression y, unsigned int mu, unsigned int lambda, int CP){
  auto M00 = fcn::sin(M_PI * lambda * x);
  auto M01 = fcn::sin(M_PI * lambda * y);
  auto M10 = fcn::sin(M_PI * mu * x);
  auto M11 = fcn::sin(M_PI * mu * y);
  return M00 * M11 + CP*M10*M01;
}

Expression chebychevt(Expression x, unsigned int order){
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
        output = 2 * x * chebychevt(x, order -1) - chebychevt(x, order - 2);
    }
    return output;
}


Expression legendret(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x;
    }
    else {
        output = (2 * order - 1)/order * x * legendret(x, order - 1) - (order - 1)/order * legendret(x, order - 2);
    }
    return output;
}
    
Expression besselt(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output=x + 1;
    }
    else{
        output = (2 * order - 1) * x * besselt(x, order - 1) + besselt(x, order - 2);
        
    }
    return output;
}

Expression laguerret(Expression x, unsigned int order){
    Expression output =0;
    if (order==0){
        output=1;
    }
    else if (order==1){
        output= 1 - x;
    }
    else{
        output = ( (2 * order  - 1 - x) * laguerret(x, order - 1) - (order - 1) * laguerret(x, order - 2))/order;
    }
    return output;
}


class tCoherentSum {
    public:
      //Takes amplitudes - 
        tCoherentSum();
        tCoherentSum(const EventType& type1, const MinuitParameterSet& mps);
        virtual ~tCoherentSum()=default;


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

    complex_t gplus(const Event& event) const{
        if (m_debug){
            INFO("Can we get time?");
            auto t = event[event.size() - 1]; //End of Event stores the time value.
            INFO("t = "<<t);

        }
        auto t = event[event.size() - 1]; //End of Event stores the time value.
        auto x = m_mps["tCoherentSum::x"]->mean(); //x mixing parameter
        auto y = m_mps["tCoherentSum::y"]->mean();

        //Mother properties
        auto D0Props = ParticleProperties::get(m_integratorAC.events().eventType().mother());
        auto D0Width = D0Props->width();
        auto D0Lifetime = D0Props->lifetime();
        auto D0Mass = D0Props->mass();

        //Rescale t -> t/2tau
        auto tp = t/(2 * D0Lifetime);

        //cosh/sinhytp
        auto _cosh = (fcn::exp(y * tp) + fcn::exp(-y * tp))/2.;
        auto _sinh =  (fcn::exp(y * tp) - fcn::exp(-y * tp))/2.;

        //cos/sin xtp

        auto _cos = fcn::cos(x * tp);
        auto _sin = fcn::sin(y * tp);

        auto val = _cos * _cosh + complex_t(0, 1) * _sin * _sinh;
        return val();
    }

    complex_t gminus(const Event& event) const{
        auto t = event[event.size() - 1]; //End of Event stores the time value.
        auto x = m_mps["tCoherentSum::x"]->mean(); //x mixing parameter
        auto y = m_mps["tCoherentSum::y"]->mean();

        //Mother properties
        auto D0Props = ParticleProperties::get(m_integratorAC.events().eventType().mother());
        auto D0Width = D0Props->width();
        auto D0Lifetime = D0Props->lifetime();
        auto D0Mass = D0Props->mass();

        //Rescale t -> t/2tau
        auto tp = t/(2 * D0Lifetime);

        //cosh/sinhytp
        auto _cosh = (fcn::exp(y * tp) + fcn::exp(-y * tp))/2.;
        auto _sinh =  (fcn::exp(y * tp) - fcn::exp(-y * tp))/2.;

        //cos/sin xtp

        auto _cos = fcn::cos(x * tp);
        auto _sin = fcn::sin(y * tp);

        auto val = _cos * _sinh + complex_t(0, 1) * _sin * _cosh;
        return val();
    }

    complex_t qp()const{
        auto absqp = m_mps["tCoherentSum::absqp"]->mean();
        auto phiqp =m_mps["tCoherentSum::phiqp"]->mean();
        complex_t val =  absqp * fcn::exp(complex_t(0, 1) * phiqp)();
        return val;


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
                  sum_i = sum_i + Cpij * CPSinPolyt(X(), Y(), i, j, 1) + Cmij * CPSinPolyt(X(), Y(), i, j, -1);
                }
                else if (m_polyType=="CP_legendre"){
                   double Cpij = getC(i,j,"P");
                   double Cmij = getC(i,j,"M");
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i + Cpij* legendret(zp, i) * legendret(zm, 2*j) + Cmij * legendret(zp, i) * legendret(zm, 2*j+1);
                   

                }
                else if (m_polyType=="antiSym_legendre"){
                   double Cij = getC(i,j);

                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * legendret(zp, i) * legendret(zm, 2*j+1);
                   

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
                  sum_i = sum_i +  Cij * chebychevt(zp, i) * chebychevt(zm, 2*j+1);
                   

                }
                else if (m_polyType=="Sym_legendre"){
                   double Cij = getC(i,j);
                   auto zp = 0.5 * (X() + Y());
                   auto zm = 0.5 * (X() - Y());
                  sum_i = sum_i +  Cij * legendret(zp, i) * legendret(zm, 2*j);
                   

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
                sum_i = sum_i + Cij * chebychevt(X(), i) * chebychevt(Y(), j); 
                }
                else if (m_polyType=="legendre"){
                    sum_i = sum_i + Cij * legendret(X(), i) * legendret(Y(), j);
                }
                else if (m_polyType=="laguerre"){
                    sum_i = sum_i + Cij * laguerret(X(), i) * laguerret(Y(), j);
                }
                else if (m_polyType=="bessel"){
                    sum_i = sum_i + Cij * besselt(X(), i) * besselt(Y(), j);
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

   
    protected:
        double  m_norm  =    {0};
        Bilinears m_normalisationsAC;

        Bilinears m_normalisationsAA;

        Bilinears m_normalisationsCC;

        MinuitParameterSet m_mps;

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



        bool m_sameTag;
        bool m_flat;

  };
}
#endif

