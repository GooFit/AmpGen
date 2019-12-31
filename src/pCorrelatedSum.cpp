#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include <limits>
using namespace AmpGen;

pCorrelatedSum::pCorrelatedSum() = default; 

pCorrelatedSum::pCorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps):
  m_A(type1, mps),
  m_B(type2.conj(true), mps),
  m_C(type1.conj(true), mps),
  m_D(type2, mps),
  m_debug(NamedParameter<bool>("pCorrelatedSum::debug", false, "Print Debug messages for pCorrelatedSum")),
  m_debugFreq(NamedParameter<int>("pCorrelatedSum::debugFreq", 1000, "Print Debug messages for pCorrelatedSum")),
  m_pdebug(NamedParameter<bool>("pCorrelatedSum::pdebug", false, "Print Debug messages for pCorrelatedSum")),
  m_pNorm(NamedParameter<bool>("pCorrelatedSum::pNorm", true, "Perform the slower normalisation for polynomial")),
  m_updateNorms(NamedParameter<bool>("pCorrelatedSum::updateNorms", true, "Update the individual norms for A,B,C,D")),
  m_order(NamedParameter<unsigned int>( "pCorrelatedSum::Order" )),
  m_polyType(NamedParameter<std::string>("pCorrelatedSum::PolyType", "simple")),
  m_mps(mps),
  m_coherentIntegral(NamedParameter<bool>("pCorrelatedSum::coherentIntegral", false, "Do the super slow method to calculate normalisation")),
  m_coherentIntegralA(NamedParameter<bool>("pCorrelatedSum::coherentIntegralA", false, "Do the super slow method to calculate normalisation")),
  m_coherentIntegralB(NamedParameter<bool>("pCorrelatedSum::coherentIntegralB", false, "Do the super slow method to calculate normalisation")),
  m_coherentIntegralC(NamedParameter<bool>("pCorrelatedSum::coherentIntegralC", false, "Do the super slow method to calculate normalisation")),
  m_coherentIntegralD(NamedParameter<bool>("pCorrelatedSum::coherentIntegralD", false, "Do the super slow method to calculate normalisation"))
{
  m_normalisationsAC.resize( m_A.matrixElements().size(), m_C.matrixElements().size() ); 
  m_normalisationsBD.resize( m_B.matrixElements().size(), m_D.matrixElements().size() ); 
  if (m_coherentIntegral){
    INFO("Using Slow option - set pCorrelatedSum::coherentIntegral false to avoid this behaviour");
    m_normalisationsAA.resize( m_A.matrixElements().size(), m_A.matrixElements().size() ); 
    m_normalisationsBB.resize( m_B.matrixElements().size(), m_B.matrixElements().size() ); 
    m_normalisationsCC.resize( m_C.matrixElements().size(), m_C.matrixElements().size() ); 
    m_normalisationsDD.resize( m_D.matrixElements().size(), m_D.matrixElements().size() ); 
  }
  if (m_coherentIntegralA) {
    INFO("Using Slow option for amplitude A");
    m_normalisationsAA.resize( m_A.matrixElements().size(), m_A.matrixElements().size() ); 
  }
  if (m_coherentIntegralB) {
      INFO("Using Slow option for amplitude B");
      m_normalisationsBB.resize( m_B.matrixElements().size(), m_B.matrixElements().size() ); 
  }
  if (m_coherentIntegralC) {
      INFO("Using Slow option for amplitude C");
      m_normalisationsCC.resize( m_C.matrixElements().size(), m_C.matrixElements().size() ); 
  }
  if (m_coherentIntegralD) {
      INFO("Using Slow option for amplitude D");
      m_normalisationsDD.resize( m_D.matrixElements().size(), m_D.matrixElements().size() ); 
  }  
}
void pCorrelatedSum::prepare(){
  if (m_debug) INFO("Preparing A");
  m_A.prepare();
  if (m_debug) INFO("Preparing B");
  m_B.prepare();
  if (m_debug) INFO("Preparing C");
  m_C.prepare();
  if (m_debug) INFO("Preparing D");
  m_D.prepare();


  if (m_debug) INFO("Preparing pCorrelatedSum");
  std::vector<size_t> changedPdfIndicesA;
  std::vector<size_t> changedPdfIndicesB;
  std::vector<size_t> changedPdfIndicesC;
  std::vector<size_t> changedPdfIndicesD;
  if (m_debug) INFO("Get Matrix Elements");
  auto matrixElementsA = m_A.matrixElements();
  auto matrixElementsB = m_B.matrixElements();
  auto matrixElementsC = m_C.matrixElements();
  auto matrixElementsD = m_D.matrixElements();
  ProfileClock clockEval; 
  if (m_debug) INFO("Call Transfer Parameters");
  m_A.transferParameters();
  m_B.transferParameters();
  m_C.transferParameters();
  m_D.transferParameters();
  if (m_debug) INFO("Get changedPdfIndiciesA for "<<matrixElementsA.size()<<" Elements");
  for ( size_t i = 0; i < matrixElementsA.size(); ++i ) {
    auto& amp = matrixElementsA[i].amp;
    if (m_debug) INFO("Prepare matrix element "<<i<<" with initial value "<<matrixElementsA[i].addressData);
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    ProfileClock clockThisElement;
    if ( m_events1 != nullptr ) {
      if ( matrixElementsA[i].addressData == 999 ) {
        if (m_debug) INFO("Set addressData for the matrixElements using data from m_events1 "<<m_events1->size());
        matrixElementsA[i].addressData = m_events1->registerExpression( amp );
      }
      if (m_debug) INFO("Update Cache for matrixElements");
      m_events1->updateCache( amp, matrixElementsA[i].addressData );
    } 
    clockThisElement.stop();
    changedPdfIndicesA.push_back(i);
    amp.resetExternals();
  }

  if (m_debug) INFO("Get changedPdfIndiciesB for "<<matrixElementsB.size()<<" Elements");
  for ( size_t i = 0; i < matrixElementsB.size(); ++i ) {
    auto& amp = matrixElementsB[i].amp;
    if (m_debug) INFO("Prepare matrix element "<<i<<" with initial value "<<matrixElementsB[i].addressData);
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    ProfileClock clockThisElement;
    if ( m_events2 != nullptr ) {
      if ( matrixElementsB[i].addressData == 999 ) {
        if (m_debug) INFO("Set addressData for the matrixElements using data from m_events2 "<<m_events2->size());
        matrixElementsB[i].addressData = m_events2->registerExpression( amp );
      }
      if (m_debug) INFO("Update Cache for matrixElements");
      m_events2->updateCache( amp, matrixElementsB[i].addressData );
    } 
    clockThisElement.stop();
    changedPdfIndicesB.push_back(i);
    amp.resetExternals();
  }

  if (m_debug) INFO("Get changedPdfIndiciesC");
  for ( size_t i = 0; i < matrixElementsC.size(); ++i ) {
    auto& amp = matrixElementsC[i].amp;
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    ProfileClock clockThisElement;
    if ( m_events1 != nullptr ) {
      if ( matrixElementsC[i].addressData == 999 ) 
        matrixElementsC[i].addressData = m_events1->registerExpression( amp );
      m_events1->updateCache( amp, matrixElementsC[i].addressData );
    } 
    clockThisElement.stop();
    changedPdfIndicesC.push_back(i);
    amp.resetExternals();
  }

  if (m_debug) INFO("Get changedPdfIndiciesD");
  for ( size_t i = 0; i < matrixElementsD.size(); ++i ) {
    auto& amp = matrixElementsD[i].amp;
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    ProfileClock clockThisElement;
    if ( m_events2 != nullptr ) {
      if ( matrixElementsD[i].addressData == 999 ) 
        matrixElementsD[i].addressData = m_events2->registerExpression( amp );
      m_events2->updateCache( amp, matrixElementsD[i].addressData );
    } 
    clockThisElement.stop();
    changedPdfIndicesD.push_back(i);
    amp.resetExternals();
  }
  if (m_debug) INFO("Get Normalisation");
  updateNorms(changedPdfIndicesA, changedPdfIndicesB, changedPdfIndicesC, changedPdfIndicesD);
  if (m_pNorm){
  m_norm = slowNorm();
  }
  else {
    m_norm = norm();
  }
  m_prepareCalls++;
}

void pCorrelatedSum::updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iB,
    const std::vector<size_t>& iC, const std::vector<size_t>& iD){
  for (auto& i : iA ) m_integratorAC.prepareExpression(m_A.matrixElements()[i].amp);
  for (auto& i : iB ) m_integratorBD.prepareExpression(m_B.matrixElements()[i].amp);
  for (auto& i : iC ) m_integratorAC.prepareExpression(m_C.matrixElements()[i].amp);
  for (auto& i : iD ) m_integratorBD.prepareExpression(m_D.matrixElements()[i].amp);

  if (m_debug) INFO("UpdateNorms for component Amplitude");
  if (m_updateNorms){ 
    m_A.updateNorms(iA);
    m_B.updateNorms(iB);
    m_C.updateNorms(iC);
    m_D.updateNorms(iD);
  if (m_coherentIntegralA) for (auto& i : iA ) m_integratorAA.prepareExpression(m_A.matrixElements()[i].amp);
  if (m_coherentIntegralB) for (auto& i : iB ) m_integratorBB.prepareExpression(m_B.matrixElements()[i].amp);
  if (m_coherentIntegralC) for (auto& i : iC ) m_integratorCC.prepareExpression(m_C.matrixElements()[i].amp);
  if (m_coherentIntegralD) for (auto& i : iD ) m_integratorDD.prepareExpression(m_D.matrixElements()[i].amp);
  }
  
  std::vector<size_t> cacheIndexA;
  std::vector<size_t> cacheIndexB;
  std::vector<size_t> cacheIndexC;
  std::vector<size_t> cacheIndexD;

  if (m_debug) INFO("integratorAC has "<<m_integratorAC.events().size());
  if (m_debug) INFO("integratorBD has "<<m_integratorBD.events().size());

  if (m_debug) INFO("Beginning to get CacheIndex for amplitudes A,B,C,D - for the integrator to use");
  if (m_debug) INFO("A index");
  for (auto& m : m_A.matrixElements() ){
    cacheIndexA.push_back(m_integratorAC.events().getCacheIndex(m.amp));
  }
  for (auto& m : m_B.matrixElements() ){
    cacheIndexB.push_back(m_integratorBD.events().getCacheIndex(m.amp));
  }
  for (auto& m : m_C.matrixElements()){
    cacheIndexC.push_back(m_integratorAC.events().getCacheIndex(m.amp));
  }
  for (auto& m : m_D.matrixElements()){
    cacheIndexD.push_back(m_integratorBD.events().getCacheIndex(m.amp));
  }
  if (m_debug) INFO("Queuing the integral for nAC"); 
  for (auto i : iA ){
    for (auto j : iC){
      m_integratorAC.queueIntegral(cacheIndexA[i], cacheIndexC[j], i, j, &m_normalisationsAC, false);
    }
  }
  if (m_debug) INFO("Finishing integral for nAC");
  m_integratorAC.flush();
  m_normalisationsAC.resetCalculateFlags();
  if (m_debug) INFO("Queuing the integral for nBD"); 
  for (auto i : iB ){
    for (auto j : iD){
      m_integratorBD.queueIntegral(cacheIndexB[i], cacheIndexD[j], i, j, &m_normalisationsBD, false);
    }
  }

  if (m_debug) INFO("Finishing integral for nBD");
  m_integratorBD.flush();
  m_normalisationsBD.resetCalculateFlags();


  if (m_coherentIntegralA){
    if (m_debug) INFO("Queuing the integral for nAA"); 
    for (auto i : iA ){
      for (auto j : iA){
        m_integratorAA.queueIntegral(cacheIndexA[i], cacheIndexA[j], i, j, &m_normalisationsAA);
      }
    }
  if (m_debug) INFO("Finishing integral for nAA");
  m_integratorAA.flush();
  m_normalisationsAA.resetCalculateFlags();
  }

  if (m_coherentIntegralB){
    if (m_debug) INFO("Queuing the integral for nBB"); 
    for (auto i : iB ){
      for (auto j : iB){
        m_integratorBB.queueIntegral(cacheIndexB[i], cacheIndexB[j], i, j, &m_normalisationsBB);
      }
    }
    if (m_debug) INFO("Finishing integral for nBB");
    m_integratorBB.flush();
    m_normalisationsBB.resetCalculateFlags();
  }

  if (m_coherentIntegralC){
    if (m_debug) INFO("Queuing the integral for nCC"); 
    for (auto i : iC ){
      for (auto j : iC){
        m_integratorCC.queueIntegral(cacheIndexC[i], cacheIndexC[j], i, j, &m_normalisationsCC);
      }
    }
    if (m_debug) INFO("Finishing integral for nCC");
    m_integratorCC.flush();
    m_normalisationsCC.resetCalculateFlags();
  }

  if (m_coherentIntegralD){
    if (m_debug) INFO("Queuing the integral for nDD"); 
    for (auto i : iD ){
      for (auto j : iD){
        m_integratorDD.queueIntegral(cacheIndexD[i], cacheIndexD[j], i, j, &m_normalisationsDD);
      }
    }
    if (m_debug) INFO("Finishing integral for nDD");
    m_integratorDD.flush();
    m_normalisationsDD.resetCalculateFlags();
  }
}

real_t pCorrelatedSum::norm() const {
  if (m_debug) INFO("Getting the value for the normalised pdf");
  if (m_debug) INFO("Get Matrix Elements for A,B,C,D");
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
  complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );

  complex_t nBD = sum_amps( m_normalisationsBD, m_B.matrixElements(), m_D.matrixElements() );
  real_t nA = m_A.norm();
  real_t nB = m_B.norm();
  real_t nC = m_C.norm();
  real_t nD = m_D.norm();

  if (m_coherentIntegralA) {
      complex_t nAA = sum_amps( m_normalisationsAA, m_A.matrixElements(), m_A.matrixElements() );
      nA = nAA.real(); 
  }
  if (m_coherentIntegralB) {
      complex_t nBB = sum_amps( m_normalisationsBB, m_B.matrixElements(), m_B.matrixElements() );
      nB = nBB.real(); 
  }
  if (m_coherentIntegralC) {
      complex_t nCC = sum_amps( m_normalisationsCC, m_C.matrixElements(), m_C.matrixElements() );
      nC = nCC.real(); 
  }
  if (m_coherentIntegralD) {
      complex_t nDD = sum_amps( m_normalisationsDD, m_D.matrixElements(), m_D.matrixElements() );
      nD = nDD.real(); 
  }

  
  complex_t mix = nAC * nBD;
  real_t intTerm = -2 * mix.real();

  

  if (m_debug) INFO("interference = "<<intTerm);
  real_t N = nA * nB + nC * nD + intTerm;

  if (m_debug) INFO("Normalisation = "<<N);
  return N;
}

double pCorrelatedSum::slowNorm(){

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
 // complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );
  std::complex<long double> nBD = sum_amps( m_normalisationsBD, m_B.matrixElements(), m_D.matrixElements() );


   std::complex<long double> rAC = 0; 
   auto eventsAC = m_integratorAC.events();

    for( int i=0 ; i < eventsAC.size(); i++ ) {
      std::complex<long double> ab = m_A.getVal(eventsAC[i]) * std::conj(m_C.getVal(eventsAC[i])) * correction(eventsAC[i])/(double)eventsAC.size();

    rAC = rAC + ab;
    
    }


  


  double n_fast = m_norm;


  double norm = m_A.norm() * m_B.norm() + m_C.norm() * m_D.norm();
  auto inter = nBD * rAC;
  norm = norm - 2 * std::real(inter);
return norm;
}

void pCorrelatedSum::reset(bool resetEvents){
  if (resetEvents){
    m_events1 = nullptr;
    m_events2 = nullptr;
    m_integratorAC = Integrator<10>();
    m_integratorBD = Integrator<10>();
    m_integratorAA = Integrator<10>();
    m_integratorBB = Integrator<10>();
    m_integratorCC = Integrator<10>();
    m_integratorDD = Integrator<10>();
    m_A.reset(true);
    m_B.reset(true);
    m_C.reset(true);
    m_D.reset(true);
  }
}

void pCorrelatedSum::setEvents(EventList& list1, EventList& list2){
  m_events1 = &list1;
  m_events2 = &list2;
  m_A.setEvents(list1);
  m_B.setEvents(list2);
  m_C.setEvents(list1);
  m_D.setEvents(list2);
}

void pCorrelatedSum::setMC(EventList& list1, EventList& list2){
  m_integratorAC = Integrator<10>(&list1);
  m_integratorBD = Integrator<10>(&list2);
  m_integratorAA = Integrator<10>(&list1);
  m_integratorBB = Integrator<10>(&list2);
  m_integratorCC = Integrator<10>(&list1);
  m_integratorDD = Integrator<10>(&list2);
  m_A.setMC(list1);
  m_B.setMC(list2);
  m_C.setMC(list1);
  m_D.setMC(list2);
}

complex_t pCorrelatedSum::getVal(const Event& evt1, const Event& evt2) const {
  complex_t A = m_A.getVal(evt1);
  complex_t B = m_B.getVal(evt2);
  complex_t C = m_C.getVal(evt1);
  complex_t D = m_D.getVal(evt2);
  complex_t f = correction(evt1);
  complex_t val = A  *f * B - C * D;
  //INFO("Correction  = "<<f);
  return val;
}
complex_t pCorrelatedSum::getValNoCache(const Event& evt1, const Event& evt2) const {
  complex_t A = m_A.getValNoCache(evt1);
  complex_t B = m_B.getValNoCache(evt2);
  complex_t C = m_C.getValNoCache(evt1);
  complex_t D = m_D.getValNoCache(evt2);
  complex_t f = correction(evt1);
  complex_t val = A *f* B - C * D;
  return val;
}


void pCorrelatedSum::debugNorm()
{
  double n_slow = 0.0;
  unsigned int nEvents = m_integratorAC.events().size();
  INFO("nEvents = "<<nEvents);
  for(unsigned int j = 0; j < m_integratorAC.events().size(); ++j ) {
    auto val = getVal(m_integratorAC.events()[j],  m_integratorBD.events()[j]);
    double val2 = val.real() * val.real() + val.imag() * val.imag();
    auto val3 = m_A.getVal(m_integratorAC.events()[j]) * std::conj(m_C.getVal(m_integratorAC.events()[j]));
    auto val4 = m_B.getVal(m_integratorBD.events()[j]) * std::conj(m_D.getVal(m_integratorBD.events()[j]));
    auto val5 = correction(m_integratorAC.events()[j]);
    n_slow += (double)val2/(double)nEvents;
    double frac = 100*(double)j/(double)nEvents;
    if (j% m_debugFreq == 0 && m_debug){
      INFO("At "<<j <<" of "<<nEvents<<" events");
      INFO("ABC*D* = "<<val);
      INFO("|AB-CD|^2 = "<<val2);
      INFO("AC* = "<<val3);
      INFO("BD* = "<<val4);
      INFO("exp(i f(s+,s-)="<<val5);
    }
  }
  auto integrate = []( const EventList& events, const CoherentSum& A, const CoherentSum& B, bool debug=false) -> double
  {
   double r = 0.0; 
   int N = events.size()/1; 
    for( int i=0 ; i < N; i++ ) {
      std::complex<long double> ab = A.getVal(events[i]) * std::conj(B.getVal(events[i]));
     double rab = std::real(ab)/N;
      
      
      
//    r += A.getVal(event) * std::conj( B.getVal(event) )/(double) events.size();
//    if (i % 100000 == 0){
  if (debug){
    INFO("r="<<r);


    INFO("A = "<<A.getVal(events[i]));
    INFO("B = "<<B.getVal(events[i]));
    INFO("AB* = "<<std::real(A.getVal(events[i]) * std::conj(B.getVal(events[i])))/events.size());
    INFO(i<<"/"<<events.size());
  }
  //  }
  long double L=events.size();
    
    
     r = r + rab;

   if (std::abs(rab)<L )  r = r + rab;

    }

    return r;
  };
  double nAA = std::real( integrate( m_integratorAC.events(), m_A, m_A, m_debug ) );
  double nBB = std::real( integrate( m_integratorBD.events(), m_B, m_B, m_debug ) );
  double nCC = std::real( integrate( m_integratorAC.events(), m_C, m_C, m_debug ) );
  double nDD = std::real( integrate( m_integratorBD.events(), m_D, m_D, m_debug ) );
  auto   nAC = ( integrate( m_integratorAC.events(), m_A, m_C, m_debug ) );
  auto   nBD = ( integrate( m_integratorBD.events(), m_B, m_D, m_debug ) );
  
  INFO("integral |A|^2 = "<<nAA);
  INFO("fast |A|^2 = "<<m_A.norm());
  INFO("|B| = "<<nBB);
  INFO("fast |B|^2 = "<<m_B.norm());
  INFO("|C| = "<<nCC);
  INFO("fast |C|^2 = "<<m_C.norm());
  INFO("|D| = "<<nDD);
  INFO("fast |D|^2 = "<<m_D.norm());
  INFO("ABC*D* = "<<-2*nAC * nBD);
  double I = ( nAA * nBB + nCC * nDD - 2*std::real(nAC*nBD) );
  double n_fast = m_norm;
  double difference = std::abs(n_fast - I)/n_fast * 100.;
  INFO("My value for n_fast = "<<n_fast);
  INFO("My value for n_slow = "<<I);
  INFO("The difference = "<<difference<<"%");
}


void pCorrelatedSum::debug(const Event& evt1, const Event& evt2) const 
{
  INFO( "A[x, y] = " << getVal(evt1, evt2) << " " << getValNoCache(evt1, evt2) << " " << m_norm << " " << prob(evt1,evt2) << " " << prob_unnormalised(evt1, evt2) );
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
