#include "AmpGen/CoherentSum.h"
#include "AmpGen/FitFraction.h"

#ifndef CORRELATEDSUM
#define CORRELATEDSUM

/* 
   CorrelatedSum is the Quantum Correlated sum, designed as an input for the QcFitter program. 
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
  class CorrelatedSum {
    public:
      //Takes amplitudes - 
      CorrelatedSum();
      CorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps);
      virtual ~CorrelatedSum()=default;
      real_t operator()( const Event& event1, const Event& event2) const { return prob(event1, event2); }
      real_t prob(const Event& event1, const Event& event2) const {
            double P = std::norm(getVal(event1, event2))/m_norm;  
            double A2 = m_A.prob(event1);
            double B2 = m_B.prob(event2);
            double C2 = m_C.prob(event1);
            double D2 = m_D.prob(event2);
            if (m_debug){
              INFO("|AB-CD|^2 = "<<P);
              INFO("|A|^2 = "<<A2);
              INFO("|B|^2 = "<<B2);
              INFO("|C|^2 = "<<C2);
              INFO("|D|^2 = "<<D2);
            }

        
        //return std::norm(getVal(event1, event2))/m_norm;
        return P;
        }
      real_t prob_unnormalised(const Event& event1, const Event& event2) const {return std::norm(getVal(event1, event2));}
      void prepare();
      void reset(bool resetEvents);

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
        std::vector<complex_t> vals = {A,B,C,D,ABCD};
        return vals;
      }
      //real_t norm(const Bilinears& norms) const; 
      double m_inter = 0;
      void updateNorm(){
        m_norm = norm();
      }
    
real_t LL(){
    real_t _LL = 0;
    real_t n = norm();
    INFO("n = "<<n);
    #pragma omp parallel for reduction( +: _LL )
    for (size_t i=0; i < m_events1->size(); i++){
        _LL += log(std::norm(getVal((*m_events1)[i], (*m_events2)[i]))/n);
    }
    return -2 * _LL;
}


    protected:
      double  m_norm  =    {0};
      Bilinears m_normalisationsAC;
      Bilinears m_normalisationsBD;

      Bilinears m_normalisationsAA;
      Bilinears m_normalisationsBB;
      Bilinears m_normalisationsCC;
      Bilinears m_normalisationsDD;

      CoherentSum  m_A;
      CoherentSum  m_B;
      CoherentSum  m_C;
      CoherentSum  m_D;
      EventList* m_events1 = nullptr;
      EventList* m_events2 = nullptr;
      Integrator<10> m_integratorAC;
      Integrator<10> m_integratorBD;



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
      size_t m_prepareCalls = 0;
  };
}
#endif
