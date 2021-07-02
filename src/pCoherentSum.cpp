#include "AmpGen/pCoherentSum.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include <limits>


#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace AmpGen;

pCoherentSum::pCoherentSum() = default; 

pCoherentSum::pCoherentSum(const EventType& type1, const MinuitParameterSet& mps,
int gammaSign, bool useXY, bool BConj):
  m_debug(NamedParameter<bool>("pCoherentSum::debug", false, "Print Debug messages for pCoherentSum")),
  
m_type1(type1), 


  m_slowNorm(NamedParameter<bool>("pCoherentSum::slowNorm", true, "Perform the slower normalisation for polynomial")),


  m_orderMag(NamedParameter<unsigned int>( "pCoherentSum::OrderMag" )),

  m_polyTypeMag(NamedParameter<std::string>("pCoherentSum::PolyTypeMag", "simple")),
  m_mps(mps),
  m_pc1(mps),
  m_pcMC1(mps),
  m_gammaSign(gammaSign),
  m_useXY(useXY),
  m_BConj(BConj),



  m_flat(NamedParameter<bool>("pCoherentSum::flat", false, "Force Amplitude=1"))





  
{

//  if (m_BConj) m_type1 = m_type1.conj(true);

  m_A = new CoherentSum(m_type1, m_mps),

  m_C = new CoherentSum(m_type1.conj(NamedParameter<bool>("pCoherentSum::ConjHead", true, "Only Conjugate the head of the EventType")), m_mps),




  INFO("Type 1 "<<type1);

}
void pCoherentSum::prepare(){

  m_A->prepare();

  m_C->prepare();

}

real_t pCoherentSum::norm() const {
  if (m_debug) INFO("Getting the value for the normalised pdf");
  complex_t sumFactor = getSumFactor(); 
  m_A->transferParameters();
  m_C->transferParameters();
  complex_t z(0,0);
  real_t w=0;
  real_t n=0;
 #pragma omp parallel for reduction( +: n )
 for (size_t i=0; i < m_sim1->size(); i++){
    //real_t f = m_pcMC1.getValCache((m_sim1)[i].address());
//    real_t f = m_pcMC1.calcCorrL((*m_sim1)[i]);
    n+=std::norm(getVal((*m_sim1)[i])); 
/*
    if (m_BConj){
      
      //z += (m_C.getValNoCache((m_sim1)[i]) * std::conj(m_A.getValNoCache((m_sim1)[i])) * exp(complex_t(0,f)));
      w += std::real((m_C.getValNoCache((*m_sim1)[i]) * std::conj(m_A.getValNoCache((*m_sim1)[i])) * exp(complex_t(0,f))) * std::conj(sumFactor))/(real_t)m_sim1->size() ;

    }
    else{
      //z += (m_A.getValNoCache((m_sim1)[i]) * std::conj(m_C.getValNoCache((m_sim1)[i])) * exp(complex_t(0,f)));
      w += std::real((m_A.getValNoCache((*m_sim1)[i]) * std::conj(m_C.getValNoCache((->m_sim1)[i])) * exp(complex_t(0,f))) * std::conj(sumFactor))/(real_t)m_sim1->size() ;
    }
    */

 }
 return n/(real_t)m_sim1->size();
  //z = z/(real_t)m_sim1.size(); 
/*
if (m_BConj){
 //return m_A.norm() + std::norm(sumFactor) * m_C.norm() + 2 * std::real(z * std::conj(sumFactor));
 return m_A.norm() + std::norm(sumFactor) * m_C.norm() + 2 * w;
}
else {
 //return m_C.norm() + std::norm(sumFactor) * m_A.norm() + 2 * std::real(z * std::conj(sumFactor));
 return m_C.norm() + std::norm(sumFactor) * m_A.norm() + 2 * w;
}
*/

}





void pCoherentSum::reset(bool resetEvents){
  if (resetEvents){
    /*
    m_events1 = nullptr;
    m_events2 = nullptr;
    m_sim1 = nullptr;
    m_sim2 = nullptr;
    */
    m_A->reset(true);

    m_C->reset(true);

  }
}

void pCoherentSum::setEvents(EventList& list1){
  m_events1 = &list1;

  m_A->setEvents(list1);

  m_C->setEvents(list1);

 // m_pc1.setEvents(list1);
 // m_pc1.prepareCache();



}

void pCoherentSum::setMC(EventList& list1){
  m_sim1 = &list1;

  m_A->setMC(list1);

  m_C->setMC(list1);

  //m_pcMC1.setEvents(list1);
  //m_pcMC1.prepareCache();
}
void pCoherentSum::setMC1(EventList& list1){
  m_sim1 = &list1;
  m_A->setMC(list1);

  m_C->setMC(list1);

}

complex_t pCoherentSum::getVal(const Event& evt1) const {
  if (m_flat) return 1;
  //auto f = m_pc1.getValCache(evt1)/2;
  auto f = m_pc1.calcCorrL(evt1)/2;
  m_A->transferParameters();


 // return m_A.getVal(evt1) * m_B.getVal(evt2) * exp( std::complex<double>(0,1) ) + getSumFactor() *  m_C.getVal(evt1) * m_D.getVal(evt2) * exp(-std::complex<double>(0,1) );
 if (m_BConj){
  return m_C->getVal(evt1) * exp(complex_t(0, f)) + getSumFactor() *  m_A->getVal(evt1) * exp(complex_t(0, -f));
 }
  else{
  return m_A->getVal(evt1) * exp(complex_t(0, f)) + getSumFactor() *  m_C->getVal(evt1) * exp(complex_t(0, -f));
 }

}
complex_t pCoherentSum::getValNoCache(const Event& evt1) const {
  complex_t A = m_A->getValNoCache(evt1);

  complex_t C = m_C->getValNoCache(evt1);

  complex_t f = correction(evt1);

  auto i = Constant(0,1);

  //auto sumFactor = getSumFactor();
  double sumFactor = -1;
  complex_t val = A  *exp(i()*f/2.)  + sumFactor * C  *exp(-i()*f/2.);
  //complex_t val = A *exp(i()*f)* B - C * D;
  if (m_flat){
    val = 1;
  }
  return val;
}


void pCoherentSum::debugNorm()
{
auto i = Constant(0,1);
  auto n_slow = 0.0;




}


void pCoherentSum::debug(const Event& evt1) const 
{
  INFO( "A[x, y] = " << getVal(evt1) << " " << getValNoCache(evt1) << " " << m_norm << " " << prob(evt1) << " " << prob_unnormalised(evt1) );
  INFO("|AB-CD|^2 = "<<prob(evt1));
  INFO("|A|^2 = "<<m_A->prob(evt1));

  INFO("|C|^2 = "<<m_C->prob(evt1));

  INFO("Norm = "<<m_norm);
}

std::vector<std::vector<FitFraction> > pCoherentSum::fitFractions(const LinearErrorPropagator& linProp){
    outputFractionsA = m_A->fitFractions(linProp);

    outputFractionsC = m_C->fitFractions(linProp);

    std::vector<std::vector<FitFraction> > outputFractions;   
    outputFractions.push_back(outputFractionsA);

    outputFractions.push_back(outputFractionsC);

    return outputFractions;
}
