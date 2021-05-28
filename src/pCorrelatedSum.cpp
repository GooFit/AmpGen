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
  m_A = CoherentSum(type1, mps);
  m_B = CoherentSum(type2.conj(NamedParameter<bool>("pCorrelatedSum::ConjHead", true, "Only Conjugate the head of the EventType")), mps);
  m_C = CoherentSum(type1.conj(NamedParameter<bool>("pCorrelatedSum::ConjHead", true, "Only Conjugate the head of the EventType")), mps);
  m_D = CoherentSum(type2, mps);


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
//  INFO("preparing");
//  prepare();

}
void pCorrelatedSum::prepare(){
  if (m_debug) INFO("Preparing A");
  m_A.transferParameters();
  m_A.prepare();
  if (m_debug) INFO("Preparing B");

  m_B.transferParameters();
  m_B.prepare();
  if (m_debug) INFO("Preparing C");

  m_C.transferParameters();
  m_C.prepare();
  if (m_debug) INFO("Preparing D");

  m_D.transferParameters();
  m_D.prepare();
  prepareACBD();
  //updateNorm();

}

real_t pCorrelatedSum::norm() const {
  if (m_debug){
 INFO("Getting the value for the normalised pdf");
// INFO("AMC(0) = "<<m_A->getValNoCache((*m_sim1)[0]));
 INFO("normA = "<<m_A.norm());
 INFO("normB = "<<m_B.norm());
 INFO("normC = "<<m_C.norm());
 INFO("normD = "<<m_D.norm());
  }
  //complex_t sumFactor = getSumFactor(); 




 if (m_analyticNorm){
  complex_t z(0,0);
  complex_t z1, z2;
 for (size_t i=0; i < (*m_sim1).size(); i++){
    //real_t f = m_pcMC1.getValCache(((*m_sim1))[i].address());
    real_t f = m_pcMC1.calcCorrL((*m_sim1)[i]);
    
    //if (m_sameTag) f = (m_pcMC1.getValCache(((*m_sim1))[i].address()) - m_pcMC2.getValCache(((*m_sim2))[i].address()));
    if (m_sameTag) f = (m_pcMC1.calcCorrL((*m_sim1)[i]) - m_pcMC2.calcCorrL((*m_sim2)[i]));

//    z += (e_A.getValNoCache(((*m_sim1))[i]) * m_B.getValNoCache(((*m_sim2))[i]) * std::conj(m_C.getValNoCache(((*m_sim1))[i]) * m_D.getValNoCache(((*m_sim2))[i])) * exp(complex_t(0,f)))/(real_t)(*m_sim1).size(); ;
    if (m_sameTag){
     // z += (m_A->getVal((*m_sim1)[i]) * m_B->getVal((*m_sim2)[i]) * std::conj(m_C->getVal((*m_sim1)[i]) * m_D->getVal((*m_sim2)[i])) * exp(complex_t(0,f)))/(real_t)(*m_sim1).size(); 
      z += m_A.getVal((*m_sim1)[i]) * std::conj(m_C.getVal((*m_sim1)[i]) * exp(complex_t(0,f))) / (real_t) (*m_sim1).size();
    }
    
    else{

      z1 += m_A.getVal((*m_sim1)[i]) *std::conj(m_C.getVal((*m_sim1)[i]) * exp(complex_t(0,f))) / (real_t) (*m_sim1).size();
      z2 += m_B.getVal((*m_sim2)[i]) *  std::conj(m_D.getVal((*m_sim2)[i])) / (real_t) (*m_sim1).size();

    }
//    z += (m_A.getVal(((*m_sim1))[i]) * m_B.getVal(((*m_sim2))[i]) * std::conj(m_C.getVal(((*m_sim1))[i]) * m_D.getVal(((*m_sim2))[i])) * exp(complex_t(0,f)));

 }
 if (m_sameTag){
   return m_A.norm()*m_B.norm() + m_C.norm() * m_D.norm() - 2 * std::real(std::norm(z));
 }
 else{ 
   z = z1 * z2;
   return m_A.norm()*m_B.norm() + m_C.norm() * m_D.norm() - 2 * std::real(z);
 }

}
else {
  real_t n=0;

  
if (m_debug){
INFO("v = "<<getValNoCache( (*m_sim1)[0], (*m_sim2)[0] ));
  }
  real_t acR, acI, bdR, bdI;
  for (size_t i=0; i<(*m_sim1).size(); i++){
    n+=std::norm(getVal((*m_sim1)[i], (*m_sim2)[i]))/(real_t) (*m_sim1).size();
  
   // complex_t ac = m_ACstMC[i] * exp(complex_t(0, m_pc1.calcCorrL((*m_sim1)[i])))/(real_t) (*m_sim1).size() ;
    //complex_t bd = m_BDstMC[i]/(real_t) (*m_sim1).size();
   // acR += std::real(ac);
   // acI += std::imag(ac);
   // bdR += std::real(bd);
   // bdI += std::imag(bd);
  }

//  complex_t AC(acR, acI);
//  complex_t BD(bdR, bdI);
//  complex_t z = AC * BD;
//  if (m_sameTag){
//    z = AC * std::conj(AC);
//  }
 // n = m_A.norm() * m_B.norm() + m_C.norm() * m_D.norm() - 2 * std::real(z);
  if (m_debug) INFO("n = "<<n);

 
  return n;
}
}





void pCorrelatedSum::reset(bool resetEvents){
  if (resetEvents){
    /*
    (*m_events1) = nullptr;
    (*m_events2) = nullptr;
    (*m_sim1) = nullptr;
    (*m_sim2) = nullptr;
    */
    m_A.reset(true);
    m_B.reset(true);
    m_C.reset(true);
    m_D.reset(true);
  }
}

void pCorrelatedSum::setEvents(EventList& list1, EventList& list2){
  m_events1= &list1;
  m_events2 = &list2;
  m_A.setEvents(list1);
  m_B.setEvents(list2);
  m_C.setEvents(list1);
  m_D.setEvents(list2);
  //m_pc1.setEvents(list1);
  //m_pc1.prepareCache();


//  if (m_sameTag) {m_pc2.setEvents(list2); m_pc2.prepareCache();}
}

void pCorrelatedSum::setMC(EventList& list1, EventList& list2){
  m_sim1 = &list1;
  m_sim2 = &list2;
  m_A.setMC(list1);
  m_B.setMC(list2);
  m_C.setMC(list1);
  m_D.setMC(list2);

  
  //m_pcMC1.setEvents(list1);
  //m_pcMC1.prepareCache();
  //if (m_sameTag) {
  //    m_pcMC2.setEvents(list2);
  //    m_pcMC2.prepareCache();
 // }

}
void pCorrelatedSum::setMC1(EventList& list1){
  (*m_sim1) = list1;
  m_A.setMC(list1);

  m_C.setMC(list1);

}
void pCorrelatedSum::setMC2(EventList& list2){
  (*m_sim2) = list2;


  m_B.setMC(list2);

  m_D.setMC(list2);
}


complex_t pCorrelatedSum::getVal(const Event& evt1, const Event& evt2) const {
  if (m_flat) return 1;
  auto f = m_pc1.calcCorrL(evt1)/2;
  if (m_debug){
    INFO("A(evt) = "<<m_A.getVal(evt1));
    INFO("B(evt) = "<<m_B.getVal(evt2));
    INFO("C(evt) = "<<m_C.getVal(evt1));
    INFO("D(evt) = "<<m_D.getVal(evt2));
  }
  //auto f = m_pc1.getValCache(evt1)/2;
  //if (m_sameTag) f = (m_pc1.getValCache(evt1) - m_pc2.getValCache(evt2))/2;
  //if (m_sameTag) f = (m_pc1.getValCache(evt1) - m_pc2.getValCache(evt2))/2;
  if (m_sameTag) f = (m_pc1.calcCorrL(evt1) - m_pc2.calcCorrL(evt2))/2;
 // return m_A.getVal(evt1) * m_B.getVal(evt2) * exp( std::complex<double>(0,1) ) + getSumFactor() *  m_C.getVal(evt1) * m_D.getVal(evt2) * exp(-std::complex<double>(0,1) );
  return m_A.getVal(evt1) * m_B.getVal(evt2) * exp(complex_t(0, f/2)) + getSumFactor() *  m_C.getVal(evt1) * m_D.getVal(evt2) * exp(complex_t(0, -f/2));
}
complex_t pCorrelatedSum::getValNoCache(const Event& evt1, const Event& evt2) const {
  complex_t A = m_A.getValNoCache(evt1);
  complex_t B = m_B.getValNoCache(evt2);
  complex_t C = m_C.getValNoCache(evt1);
  complex_t D = m_D.getValNoCache(evt2);
  complex_t f = m_pc1.calcCorrL(evt1);
  if (m_sameTag){
    f -= m_pc1.calcCorrL(evt2);
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

//INFO("A0 = "<<m_A.getVal((*m_sim1)[0]));
INFO("ACst[0] = "<<m_ACstMC[0]);
  real_t nA_0 = m_A.norm();
  real_t nB_0 = m_B.norm();
  real_t nC_0 = m_C.norm();
  real_t nD_0 = m_D.norm();
  real_t nA_1 = 0;
  for (size_t i=0;i<(*m_sim1).size();i++){
    nA_1 += std::norm(m_A.getValNoCache((*m_sim1)[i]));
  }
  nA_1 = nA_1/(real_t)(*m_sim1).size();
  INFO("nA_0 = "<<nA_0);
  INFO("nB_0 = "<<nB_0);
  INFO("nC_0 = "<<nC_0);
  INFO("nD_0 = "<<nD_0);
  INFO("nA_1 = "<<nA_1);

  real_t n0 = norm();
  real_t n1 = 0;
  complex_t in1 = 0;
  complex_t in2= 0;
  for (size_t i=0;i<(*m_sim1).size();i++){
    n1 += std::norm(getValNoCache((*m_sim1)[i], (*m_sim2)[i]));
    in1 += m_A.getValNoCache((*m_sim1)[i]) * std::conj(m_C.getValNoCache((*m_sim1)[i])) * exp(complex_t(0, m_pc1.calcCorrL((*m_sim1)[i])))/(real_t)(*m_sim1).size();
    in2 += m_B.getValNoCache((*m_sim2)[i]) * std::conj(m_D.getValNoCache((*m_sim2)[i]))/(real_t)(*m_sim1).size();

  }
  
  n1 = n1/(real_t)(*m_sim1).size();
  real_t z0 = n0 - m_A.norm() * m_B.norm()- m_C.norm() * m_D.norm();
  real_t z1 = n1 - m_A.norm() * m_B.norm()- m_C.norm() * m_D.norm();

  real_t z2 =-2 * std::real( in1 * in2);
  

  INFO("n0 = "<<n0);
  INFO("n1 = "<<n1);

  INFO("z1 = "<<z0);
  INFO("z1 = "<<z1);
  INFO("z2 = "<<z2);

  updateNorm();
  real_t pA =0;
  real_t p=0;
  for (size_t i=0;i<(*m_events1).size();i++){
    pA+=std::norm(m_A.getValNoCache((*m_events1)[i]))/m_A.norm();
    p+=prob((*m_events1)[i], (*m_events2)[i]);
  }
  INFO("sum prob(A) = "<<pA);

  

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
  INFO("|A|^2 = "<<m_A.prob(evt1));
  INFO("|B|^2 = "<<m_B.prob(evt2));
  INFO("|C|^2 = "<<m_C.prob(evt1));
  INFO("|D|^2 = "<<m_D.prob(evt2));
  INFO("Norm = "<<m_norm);
}

std::vector<std::vector<FitFraction> > pCorrelatedSum::fitFractions(const LinearErrorPropagator& linProp){
    outputFractionsA = m_A.fitFractions(linProp);
    outputFractionsB = m_B.fitFractions(linProp);
    outputFractionsC = m_C.fitFractions(linProp);
    outputFractionsD = m_D.fitFractions(linProp);
    std::vector<std::vector<FitFraction> > outputFractions;   
    outputFractions.push_back(outputFractionsA);
    outputFractions.push_back(outputFractionsB);
    outputFractions.push_back(outputFractionsC);
    outputFractions.push_back(outputFractionsD);
    return outputFractions;
}
