#ifndef AMPGEN_gCOHERENTSUM_H
#define AMPGEN_gCOHERENTSUM_H

#include <memory.h>
#include <stddef.h>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/CoherentSum.h"
#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Integrator2.h"
#include "AmpGen/Types.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
//#include "AmpGen/functional/pdf.h"

namespace AmpGen
{

class gCoherentSum // : public functional::pdf_base<CoherentSum>
  {
  public:
    gCoherentSum();
    gCoherentSum( const EventType& type, const AmpGen::MinuitParameterSet& mps);
    

    virtual ~gCoherentSum() = default; 
    double norm() const {


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
      
      auto nAC = sum_amps(m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements()); 
      auto inter = nAC * complex_t(m_x->mean(), -m_y->mean());
      //return m_A.norm() + m_C.norm() + 2 * inter.real();
      return m_A.norm() +  (pow(m_x->mean(), 2) + pow(m_y->mean(), 2)) * m_C.norm() + 2 * inter.real();
      //return m_C.norm() +  (pow(m_x->mean(), 2) + pow(m_y->mean(), 2)) * m_C.norm() + 2 * inter.real();
      //return m_A.norm() * (pow(m_x->mean(), 2) + pow(m_y->mean(), 2));// + m_C.norm();// + 2 * inter.real();

    }
    complex_t getVal(const Event& evt) const {
      bool flat = NamedParameter<bool>("gCoherentSum::flat", false);
      if (flat) return 1;
      return m_A.getVal(evt) + complex_t(m_x->mean(), m_y->mean()) * m_C.getVal(evt);
      //return m_A.getVal(evt) * complex_t(m_x->mean(), m_y->mean());// * m_C.getVal(evt);
//      return m_A.getVal(evt);// + complex_t(m_x->mean(), m_y->mean()) * m_C.getVal(evt);
    }

    void setEvents(EventList& list){
      m_A.setEvents(list); 
      m_C.setEvents(list); 
      m_events = &list;
    }

    void setMC(EventList& sim){
      m_A.setMC(sim);
      m_C.setMC(sim);
      m_integrator = integrator(&sim);
    }

    real_t prob(const Event& evt) const { return std::norm(getVal(evt))/m_norm; }

    real_t prob_unnormalised( const Event& evt ) const { return std::norm(getVal(evt)); }


    size_t size() const { return m_A.matrixElements().size() + m_C.matrixElements().size(); }



    real_t operator()( const Event& evt )        const { return std::norm(getVal(evt))/m_norm; }
  void reset( bool resetEvents )
  {

  for ( auto& mE : m_A.matrixElements() ) mE.addressData = 999;
  for ( auto& mE : m_C.matrixElements() ) mE.addressData = 999;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integrator = integrator();
  }
}



    void prepare(){
      m_A.prepare(); 
      m_C.prepare(); 
      m_normalisationsAC.resize( m_A.matrixElements().size(), m_C.matrixElements().size() ); 
      m_norm = norm();
    }
  protected:
    double m_norm;

    typedef Integrator<10> integrator;
    CoherentSum m_A;
    CoherentSum m_C;
    //Real and imaginary part for complex sumFactor
    MinuitParameter* m_x;
    MinuitParameter* m_y;
    //AmpGen::MinuitParameterSet m_mps;
    EventType m_type;
    EventType m_typeConj;

    EventList * m_events;
    integrator m_integrator;
    Bilinears m_normalisationsAC;
    
  };
}
#endif