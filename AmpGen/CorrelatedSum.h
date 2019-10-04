#include "AmpGen/CoherentSum.h"
#ifndef CORRELATEDSUM
#define CORRELATEDSUM

using namespace std;
using namespace AmpGen;

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

//template <typename PDF>
class CorrelatedSum{
    public:
        //Takes amplitudes - 
        CorrelatedSum();
        CorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps);
        virtual ~CorrelatedSum()=default;
        real_t operator()( const Event& event1, const Event& event2) const { return prob(event1, event2); }
        real_t prob(const Event& event1, const Event& event2) const {return std::norm(getVal(event1, event2)/m_norm);}
        real_t prob_unnormalised(const Event& event1, const Event& event2) const {return std::norm(getVal(event1, event2));}
        void prepare();
        void reset(bool resetEvents);

        void setEvents(EventList& list1, EventList& list2);
        void setMC(EventList& list1, EventList& list2);

        void updateNorms(const std::vector<unsigned int>& iA, const std::vector<unsigned int>& iB,
                         const std::vector<unsigned int>& iC, const std::vector<unsigned int>& iD);

        void debugNorm();

        complex_t getVal(const Event& event1, const Event& event2) const;
        complex_t getValNoCache(const Event& event1, const Event& event2) const;
        complex_t getValNoCache(const Event& event1, const Event& event2, const size_t& offset) const;
        real_t norm()  const;
        //real_t norm(const Bilinears& norms) const; 
        double m_inter = 0;
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


        Integrator<10> m_integratorAA;
        Integrator<10> m_integratorBB;
        Integrator<10> m_integratorCC;
        Integrator<10> m_integratorDD;

        bool m_useCoherent = false;
        bool m_coherentIntegral = false;

        bool m_debug = false;
        size_t m_prepareCalls = 0;
};

#endif
