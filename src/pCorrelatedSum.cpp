#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include <limits>


#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace AmpGen;

pCorrelatedSum::pCorrelatedSum() = default; 

pCorrelatedSum::pCorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps, std::string SFType):
  m_debug(NamedParameter<bool>("pCorrelatedSum::debug", false, "Print Debug messages for pCorrelatedSum")),
  
m_type1(type1), 
m_type2(type2), 

  m_slowNorm(NamedParameter<bool>("pCorrelatedSum::slowNorm", true, "Perform the slower normalisation for polynomial")),


  m_orderMag(NamedParameter<unsigned int>( "pCorrelatedSum::OrderMag" )),

  m_polyTypeMag(NamedParameter<std::string>("pCorrelatedSum::PolyTypeMag", "simple")),
  m_mps(mps),
  m_pc1(mps),
  m_pcMC1(mps),
  m_pc2(mps),
  m_pcMC2(mps),

  m_flat(NamedParameter<bool>("pCorrelatedSum::flat", false, "Force Amplitude=1")),
  m_analyticNorm(NamedParameter<bool>("pCorrelatedSum::analyticNorm", false)),





  m_SFType(SFType)
{
  m_A = new CoherentSum(type1, mps);
  m_B = new CoherentSum(type2.conj(NamedParameter<bool>("pCorrelatedSum::ConjHead", true, "Only Conjugate the head of the EventType")), mps);
  m_C = new CoherentSum(type1.conj(NamedParameter<bool>("pCorrelatedSum::ConjHead", true, "Only Conjugate the head of the EventType")), mps);
  m_D = new CoherentSum(type2, mps);


  INFO("Type 1 "<<type1);
  INFO("Type 2 "<<type2);
  if (type1==type2){
    m_sameTag = true;
  }
  else{
    m_sameTag = false;
  }
  if (m_sameTag){
    INFO("Using same tag for AC and BD");
  }
}
void pCorrelatedSum::prepare(){
  if (m_debug) INFO("Preparing A");
  m_A->transferParameters();
  m_A->prepare();
  if (m_debug) INFO("Preparing B");

  m_B->transferParameters();
  m_B->prepare();
  if (m_debug) INFO("Preparing C");

  m_C->transferParameters();
  m_C->prepare();
  if (m_debug) INFO("Preparing D");

  m_D->transferParameters();
  m_D->prepare();
  updateNorm();

}

real_t pCorrelatedSum::norm() const {
  if (m_debug){
 INFO("Getting the value for the normalised pdf");
 INFO("AMC(0) = "<<m_A->getValNoCache((*m_sim1)[0]));
 INFO("normA = "<<m_A->norm());
 INFO("normB = "<<m_B->norm());
 INFO("normC = "<<m_C->norm());
 INFO("normD = "<<m_D->norm());
  }
  //complex_t sumFactor = getSumFactor(); 




 if (m_analyticNorm){
  complex_t z(0,0);
 for (size_t i=0; i < m_sim1->size(); i++){
    //real_t f = m_pcMC1.getValCache((m_sim1)[i].address());
    real_t f = m_pcMC1.calcCorrL((*m_sim1)[i]);
    //if (m_sameTag) f = (m_pcMC1.getValCache((m_sim1)[i].address()) - m_pcMC2.getValCache((m_sim2)[i].address()));
    if (m_sameTag) f = (m_pcMC1.calcCorrL((*m_sim1)[i]) - m_pcMC2.calcCorrL((*m_sim2)[i]));

//    z += (e_A.getValNoCache((m_sim1)[i]) * m_B.getValNoCache((m_sim2)[i]) * std::conj(m_C.getValNoCache((m_sim1)[i]) * m_D.getValNoCache((m_sim2)[i])) * exp(complex_t(0,f)))/(real_t)m_sim1.size(); ;
    if (m_sameTag){
     // z += (m_A->getVal((*m_sim1)[i]) * m_B->getVal((*m_sim2)[i]) * std::conj(m_C->getVal((*m_sim1)[i]) * m_D->getVal((*m_sim2)[i])) * exp(complex_t(0,f)))/(real_t)m_sim1->size(); 
      z += m_A->getVal((*m_sim1)[i]) * std::conj(m_C->getVal((*m_sim1)[i]) * exp(complex_t(0,f)));
    }
    
    else{

      z += (m_A->getVal((*m_sim1)[i]) * m_B->getVal((*m_sim2)[i]) * std::conj(m_C->getVal((*m_sim1)[i]) * m_D->getVal((*m_sim2)[i])) * exp(complex_t(0,f)));
    }
//    z += (m_A.getVal((m_sim1)[i]) * m_B.getVal((m_sim2)[i]) * std::conj(m_C.getVal((m_sim1)[i]) * m_D.getVal((m_sim2)[i])) * exp(complex_t(0,f)));

 }
 if (m_sameTag){
   return m_A->norm()*m_B->norm() + m_C->norm() * m_D->norm() - 2 * std::real(std::norm(z)/m_sim1->size());
 }
 else{ 
   return m_A->norm()*m_B->norm() + m_C->norm() * m_D->norm() - 2 * std::real(z/(real_t)m_sim1->size());
 }

}
else {
  real_t n=0;
  for (size_t i=0; i<m_sim1->size(); i++){
    n+=std::norm(getVal((*m_sim1)[i], (*m_sim2)[i]))/(real_t) m_sim1->size();
  }
  return n;
}
}




void pCorrelatedSum::reset(bool resetEvents){
  if (resetEvents){
    /*
    m_events1 = nullptr;
    m_events2 = nullptr;
    m_sim1 = nullptr;
    m_sim2 = nullptr;
    */
    m_A->reset(true);
    m_B->reset(true);
    m_C->reset(true);
    m_D->reset(true);
  }
}

void pCorrelatedSum::setEvents(EventList& list1, EventList& list2){
  m_events1 = &list1;
  m_events2 = &list2;
  m_A->setEvents(list1);
  m_B->setEvents(list2);
  m_C->setEvents(list1);
  m_D->setEvents(list2);
  //m_pc1.setEvents(list1);
  //m_pc1.prepareCache();


//  if (m_sameTag) {m_pc2.setEvents(list2); m_pc2.prepareCache();}
}

void pCorrelatedSum::setMC(EventList& list1, EventList& list2){
  m_sim1 = &list1;
  m_sim2 = &list2;
  m_A->setMC(list1);
  m_B->setMC(list2);
  m_C->setMC(list1);
  m_D->setMC(list2);
  //m_pcMC1.setEvents(list1);
  //m_pcMC1.prepareCache();
  //if (m_sameTag) {
  //    m_pcMC2.setEvents(list2);
  //    m_pcMC2.prepareCache();
 // }
}
void pCorrelatedSum::setMC1(EventList& list1){
  m_sim1 = &list1;
  m_A->setMC(list1);

  m_C->setMC(list1);

}
void pCorrelatedSum::setMC2(EventList& list2){
  m_sim2 = &list2;


  m_B->setMC(list2);

  m_D->setMC(list2);
}


complex_t pCorrelatedSum::getVal(const Event& evt1, const Event& evt2) const {
  if (m_flat) return 1;
  auto f = m_pc1.calcCorrL(evt1)/2;
  if (m_debug){
    INFO("A(evt) = "<<m_A->getVal(evt1));
    INFO("B(evt) = "<<m_B->getVal(evt2));
    INFO("C(evt) = "<<m_C->getVal(evt1));
    INFO("D(evt) = "<<m_D->getVal(evt2));
  }
  //auto f = m_pc1.getValCache(evt1)/2;
  //if (m_sameTag) f = (m_pc1.getValCache(evt1) - m_pc2.getValCache(evt2))/2;
  //if (m_sameTag) f = (m_pc1.getValCache(evt1) - m_pc2.getValCache(evt2))/2;
  if (m_sameTag) f = (m_pc1.calcCorrL(evt1) - m_pc2.calcCorrL(evt2))/2;
 // return m_A.getVal(evt1) * m_B.getVal(evt2) * exp( std::complex<double>(0,1) ) + getSumFactor() *  m_C.getVal(evt1) * m_D.getVal(evt2) * exp(-std::complex<double>(0,1) );
  return m_A->getVal(evt1) * m_B->getVal(evt2) * exp(complex_t(0, f)) + getSumFactor() *  m_C->getVal(evt1) * m_D->getVal(evt2) * exp(complex_t(0, -f));
}
complex_t pCorrelatedSum::getValNoCache(const Event& evt1, const Event& evt2) const {
  complex_t A = m_A->getValNoCache(evt1);
  complex_t B = m_B->getValNoCache(evt2);
  complex_t C = m_C->getValNoCache(evt1);
  complex_t D = m_D->getValNoCache(evt2);
  complex_t f = correction(evt1);
  if (m_sameTag){
    f -= correction(evt2);
  }
  auto i = Constant(0,1);

  //auto sumFactor = getSumFactor();
  double sumFactor = -1;
  complex_t val = A  *exp(i()*f/2.) * B + sumFactor * C * D  *exp(-i()*f/2.);
  //complex_t val = A *exp(i()*f)* B - C * D;
  if (m_flat){
    val = 1;
  }
  return val;
}


void pCorrelatedSum::debugNorm()
{
auto i = Constant(0,1);
  auto n_slow = 0.0;




}


void pCorrelatedSum::debug(const Event& evt1, const Event& evt2) const 
{
  INFO( "A[x, y] = " << getVal(evt1, evt2));
  INFO("A[x,y]_NoCache = "<<getValNoCache(evt1, evt2));
  INFO("m_Norm =  " << m_norm);
  INFO("Norm =  " << norm());
  INFO("prob = " << prob(evt1,evt2));
  INFO("prob_unorm" << prob_unnormalised(evt1, evt2) );
  INFO("|AB-CD|^2 = "<<prob(evt1, evt2));
  INFO("|A|^2 = "<<m_A->prob(evt1));
  INFO("|B|^2 = "<<m_B->prob(evt2));
  INFO("|C|^2 = "<<m_C->prob(evt1));
  INFO("|D|^2 = "<<m_D->prob(evt2));
  INFO("Norm = "<<m_norm);
}

std::vector<std::vector<FitFraction> > pCorrelatedSum::fitFractions(const LinearErrorPropagator& linProp){
    outputFractionsA = m_A->fitFractions(linProp);
    outputFractionsB = m_B->fitFractions(linProp);
    outputFractionsC = m_C->fitFractions(linProp);
    outputFractionsD = m_D->fitFractions(linProp);
    std::vector<std::vector<FitFraction> > outputFractions;   
    outputFractions.push_back(outputFractionsA);
    outputFractions.push_back(outputFractionsB);
    outputFractions.push_back(outputFractionsC);
    outputFractions.push_back(outputFractionsD);
    return outputFractions;
}
