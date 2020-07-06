#include "AmpGen/tCoherentSum.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include <limits>
#include "AmpGen/ParticleProperties.h"

#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace AmpGen;

tCoherentSum::tCoherentSum() = default; 

tCoherentSum::tCoherentSum(const EventType& type1, const MinuitParameterSet& mps, std::string SFType, int gamSign, bool useXY):
  m_A(type1, mps),

  m_C(type1.conj(NamedParameter<bool>("tCoherentSum::ConjHead", true, "Only Conjugate the head of the EventType")), mps),

  m_debug(NamedParameter<bool>("tCoherentSum::debug", false, "Print Debug messages for tCoherentSum")),
  m_debugFreq(NamedParameter<int>("tCoherentSum::debugFreq", 1000, "Print Debug messages for tCoherentSum")),
  m_pdebug(NamedParameter<bool>("tCoherentSum::pdebug", false, "Print Debug messages for tCoherentSum")),
  m_pNorm(NamedParameter<bool>("tCoherentSum::pNorm", true, "Perform the slower normalisation for polynomial")),
  m_updateNorms(NamedParameter<bool>("tCoherentSum::updateNorms", true, "Update the individual norms for A,B,C,D")),
  m_order(NamedParameter<unsigned int>( "pCorrelatedSum::Order" )),
  m_polyType(NamedParameter<std::string>("pCorrelatedSum::PolyType", "simple")),
  m_mps(mps),
  m_gamSign(gamSign),
  m_useXY(useXY),
  m_flat(NamedParameter<bool>("tCoherentSum::flat", false, "Force Amplitude=1")),
  m_coherentIntegral(NamedParameter<bool>("tCoherentSum::coherentIntegral", false, "Do the super slow method to calculate normalisation")),
  m_coherentIntegralA(NamedParameter<bool>("tCoherentSum::coherentIntegralA", false, "Do the super slow method to calculate normalisation")),
  m_fastNorm(NamedParameter<bool>("tCoherentSum::fastNorm", false)),
  m_coherentIntegralC(NamedParameter<bool>("tCoherentSum::coherentIntegralC", false, "Do the super slow method to calculate normalisation")),

  m_SFType(SFType)
{
  m_normalisationsAC.resize( m_A.matrixElements().size(), m_C.matrixElements().size() ); 

  if (m_coherentIntegral){
    INFO("Using Slow option - set tCoherentSum::coherentIntegral false to avoid this behaviour");
    m_normalisationsAA.resize( m_A.matrixElements().size(), m_A.matrixElements().size() ); 

    m_normalisationsCC.resize( m_C.matrixElements().size(), m_C.matrixElements().size() ); 

  }
  if (m_coherentIntegralA) {
    INFO("Using Slow option for amplitude A");
    m_normalisationsAA.resize( m_A.matrixElements().size(), m_A.matrixElements().size() ); 
  }

  if (m_coherentIntegralC) {
      INFO("Using Slow option for amplitude C");
      m_normalisationsCC.resize( m_C.matrixElements().size(), m_C.matrixElements().size() ); 
  }
  INFO("Type 1 "<<type1);
  
  
}
void tCoherentSum::prepare(){
  if (m_debug) INFO("Preparing A");
  m_A.prepare();
  
  if (m_debug) INFO("Preparing C");
  m_C.prepare();
 
 



  if (m_debug) INFO("Preparing tCoherentSum");
  std::vector<size_t> changedPdfIndicesA;
 
  std::vector<size_t> changedPdfIndicesC;
 
  if (m_debug) INFO("Get Matrix Elements");
  auto matrixElementsA = m_A.matrixElements();
 
  auto matrixElementsC = m_C.matrixElements();
 
  ProfileClock clockEval; 
  if (m_debug) INFO("Call Transfer Parameters");
  m_A.transferParameters();
 
  m_C.transferParameters();
 
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

 
  if (m_debug) INFO("Get Normalisation");

  
  if (m_debug) INFO("length of changedPdfIndiciesA = "<<changedPdfIndicesA.size());
  if (m_debug) INFO("length of changedPdfIndiciesC = "<<changedPdfIndicesC.size());
  if (m_debug) INFO("Is integratorAA ready "<<m_integratorAA.isReady());
  if (m_debug) INFO("Is integratorAC ready "<<m_integratorAC.isReady());
  if (m_debug) INFO("Is integratorCC ready "<<m_integratorAC.isReady());
  if (m_integratorAA.isReady() && m_integratorAC.isReady() && m_integratorCC.isReady() ) updateNorms(changedPdfIndicesA,  changedPdfIndicesC );

  
  m_Anorm = m_A.norm();
  if (m_debug) INFO("Got Normalisation for m_A");
 
 
  m_Cnorm = m_C.norm();
  if (m_debug) INFO("Got Normalisation for m_C");
 
 
  m_norm = norm();
  if (m_debug) INFO("Got Normalisation from full integral");


  m_prepareCalls++;
  if (m_debug) INFO("CoherentSum is prepared!");
}

void tCoherentSum::updateNorms(const std::vector<size_t>& iA, const std::vector<size_t>& iC){
  for (auto& i : iA ) m_integratorAC.prepareExpression(m_A.matrixElements()[i].amp);


  for (auto& i : iC ) m_integratorAC.prepareExpression(m_C.matrixElements()[i].amp);


  if (m_debug) INFO("UpdateNorms for component Amplitude");
  if (m_updateNorms){ 
    m_A.updateNorms(iA);

    m_C.updateNorms(iC);

  if (m_coherentIntegralA) for (auto& i : iA ) m_integratorAA.prepareExpression(m_A.matrixElements()[i].amp);

  if (m_coherentIntegralC) for (auto& i : iC ) m_integratorCC.prepareExpression(m_C.matrixElements()[i].amp);

  }
  
  std::vector<size_t> cacheIndexA;

  std::vector<size_t> cacheIndexC;
 // if (m_debug) INFO("Getting events from integratorAC"); 
 // auto ACevents = m_integratorAC.events();
//  if (m_debug) INFO("Got events from integratorAC"); 



  if (m_debug) INFO("integratorAC has "<<m_integratorAC.events().size());


  if (m_debug) INFO("Beginning to get CacheIndex for amplitudes A,C - for the integrator to use");
  if (m_debug) INFO("A index");
  for (auto& m : m_A.matrixElements() ){
    cacheIndexA.push_back(m_integratorAC.events().getCacheIndex(m.amp));
  }
  
  for (auto& m : m_C.matrixElements()){
    cacheIndexC.push_back(m_integratorAC.events().getCacheIndex(m.amp));
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
  

  if (m_debug) INFO("Finished Updating Norms");

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

  
}



real_t tCoherentSum::norm() const {
  if (m_debug) INFO("Getting the value for the normalised pdf");
   auto eventsAC = m_integratorAC.events();
  if (m_debug) INFO("Get sumFactor");
 
  complex_t sumFactor = getSumFactor();  
  if (m_debug) INFO("sumFactor="<<sumFactor);

  if (m_debug) INFO("Get TD params");
  if (m_debug) INFO("Get x");
   auto param_x = m_mps["tCoherentSum::x"];

  if (m_debug) INFO("x = "<<param_x->mean());
  auto param_y = m_mps["tCoherentSum::y"];
  if (m_debug) INFO("y = "<<param_y->mean());
  auto param_absqp = m_mps["tCoherentSum::absqp"];
  if (m_debug) INFO("|q/p| = "<<param_absqp->mean());
  auto param_phiqp = m_mps["tCoherentSum::phiqp"];
  if (m_debug) INFO("phi(q/p) = "<<param_phiqp->mean());

  auto evtType = eventsAC.eventType();
  if (m_debug) INFO("eventType "<<evtType);
  if (m_debug) INFO("Get mother "<<eventsAC.eventType().mother()<<" properties");

  auto D0Props = ParticleProperties::get(eventsAC.eventType().mother());
  auto D0Width = D0Props->width();
  auto D0Lifetime = D0Props->lifetime();
  if (m_debug) INFO(eventsAC.eventType().mother()<<" width = "<<D0Width);
  
  auto D0Mass = D0Props->mass();
auto h = 6.62607015e-34;
auto c = 299792458.0;
D0Mass = D0Mass * pow(c, 2)/h;

  if (m_debug) INFO(eventsAC.eventType().mother()<<" mass = "<<D0Mass);
  auto mix_x = param_x->mean();
  if (m_debug) INFO("x = "<<mix_x);
  auto mix_y = param_y->mean();
  if (m_debug) INFO("y = "<<mix_y);
  auto absqp = param_absqp->mean();
  auto phiqp = param_phiqp->mean();
  auto qp = absqp * fcn::exp(Constant(0, 1)() * phiqp)(); 

  if (m_debug) INFO("q/p = "<<qp);
//  auto a = fcn::exp(-D0Width * t - mix_y * D0Width * t/2 - i() * D0Mass * t + i * mix_x * D0Width * t/2); 
//  auto b = fcn::exp(-i() * mix_x * D0Width * t/2 + D0Width * mix_y * t);
//  auto gplus = a * (1 + b);
// auto  gminus = a * (1-b);  



//  INFO("Int |g_+(t)|^2 dt = "<<int_gp2());
  //auto int_gm2 = 1/(2/D0Lifetime + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  - 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );


//  INFO("Int |g_-(t)|^2 dt = "<<int_gp2());

//  INFO("Int g_+(t)g_-^*(t) dt = "<<int_gpm());


//  if (m_debug) INFO("Calculate TD parts of normalisation");
//  auto int_gp2 = 1/(2 * D0Width + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  + 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );
//  auto int_gp2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  + 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );

 // if (m_debug) INFO("Int |g_+(t)|^2 dt = "<<int_gp2());
  //auto int_gm2 = 1/(2 * D0Width + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  - 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );
//  auto int_gm2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  + 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
//  if (m_debug) INFO("Int |g_-(t)|^2 dt = "<<int_gp2());
//  auto int_gpm = 1/(2 * D0Width + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width) - Constant(0,1)() * mix_x * D0Width/(4 * fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)()/4);
//  auto int_gpm = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime) - Constant(0,1)() * mix_x * (1/D0Lifetime)/(4 * fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x/D0Lifetime, 2)()/4);
//  if (m_debug) INFO("Int g_+(t)g_-^*(t) dt = "<<int_gpm());

  //auto int_gp2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) + 1/(1 + fcn::pow(mix_x, 2)));
  auto int_gp2 = (1/2)  * (1/(1 - fcn::pow(mix_y, 2)) + 1/(1 + fcn::pow(mix_x, 2)));
  if (m_debug) INFO("int |g_+|^2 = "<<int_gp2());

  //auto int_gm2 = 1/(2/D0Lifetime + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  - 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );

  //auto int_gm2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * (D0Lifetime/(1 - fcn::pow(mix_y, 2)) - D0Lifetime/(1 + fcn::pow(mix_x, 2))); //  (2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) - 1/(1 + fcn::pow(mix_x, 2)));
  auto int_gm2 = (1/2) * (1/(1 - fcn::pow(mix_y, 2)) - 1/(1 + fcn::pow(mix_x, 2)));

  if (m_debug) INFO("int |g_-|^2 = "<<int_gm2());

  //auto int_gpm = D0Lifetime/(2 + mix_y) + D0Lifetime/(-2 + mix_y) - Constant(0,1)() * mix_x * D0Lifetime/(4 + fcn::pow(mix_x, 2)()/4);
  auto int_gpm = 1/(2 + mix_y) + 1/(-2 + mix_y) - Constant(0,1)() * mix_x * 1/(4 + fcn::pow(mix_x, 2)()/4);

  if (m_debug) INFO("int g_+g_-* = "<<int_gpm());

 if (m_debug) INFO("Get norms for A,C");



  real_t nA = m_A.norm();

  real_t nC = m_C.norm();


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
  if (m_debug) INFO("Get interference without corrections");
  complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );

  if (m_debug) INFO("interference = "<<nAC);

  complex_t mix = nAC * std::conj(sumFactor) * int_gpm() * qp;
  //complex_t mix = nAC * std::conj(sumFactor) * qp;
  if (m_debug) INFO("Calculate uncorrected normalisation");
  //real_t N = int_gp2().real() * nA + std::norm(sumFactor) * std::norm(qp) * int_gm2().real() * nC + 2 * std::real(int_gpm() * qp * mix);
  auto N = int_gp2() * nA + std::norm(sumFactor) * std::norm(qp) * int_gm2() * nC + 2 * mix.real();
  //auto N =  nA + std::norm(sumFactor) * nC + 2 * mix.real();
  
  if (m_fastNorm){
    return N.real();    
  }


   std::complex<double> rAC = 0; 



if (m_debug){
for (int i=0; i<eventsAC.size(); i++){
  INFO("Can I call getVal for A at "<<i<<" out of "<<eventsAC.size());
  INFO("Event = "<<eventsAC[i]);
  INFO("m_A = "<<m_A.getValNoCache(eventsAC[i]));  
  INFO("m_C = "<<m_C.getValNoCache(eventsAC[i]));  
}
}

   double N_slow = 0; 
    for( int i=0 ; i < eventsAC.size(); i++ ) {
  auto t = eventsAC[i][4 * eventsAC.eventType().size()];
  auto a = fcn::exp(- t/D0Lifetime - mix_y * t/(2*D0Lifetime) - Constant(0, 1) * D0Mass * t + Constant(0, 1) * mix_x * t/(2*D0Lifetime)); 
  auto b = fcn::exp(-1 * Constant(0,1) * mix_x * t/(D0Lifetime) + mix_y * t/D0Lifetime);
 auto gplus = a * (1 + b);
auto  gminus = a * (1-b);

    auto vi = std::norm(getValNoCache(eventsAC[i]));
    N_slow += vi;
      auto ab = std::conj(qp) * gplus() * std::conj(gminus()) * m_A.getValNoCache(eventsAC[i]) * std::conj(m_C.getValNoCache(eventsAC[i])) * exp(Constant(0,1)() * correction(eventsAC[i]))/(double)eventsAC.size();
  //    if (std::abs(ab) > 1e3){
 
      //auto abA = abs(m_A.getValNoCache(eventsAC[i]));
      auto abA = abs(m_A.getValNoCache(eventsAC[i]) * gplus());
      //auto abC = abs(m_C.getValNoCache(eventsAC[i]));
      auto abC = abs(m_C.getValNoCache(eventsAC[i]) * qp * gminus());
    
      if (abA>1e10 || abC>1e10){
  //      INFO("A or C too large?");
      }
      else{
//   if (m_debug)   INFO("ab = "<<abs(ab));
//    if (m_debug)    INFO("correction = "<<correction(eventsAC[i]));
//      if (m_debug)    INFO("rAC = "<<abs(rAC));

        rAC = rAC + ab;
    }
    }

  N_slow = N_slow/(double)eventsAC.size();

  if (m_coherentIntegralA) {
      complex_t nAA = sum_amps( m_normalisationsAA, m_A.matrixElements(), m_A.matrixElements() );
      nA = nAA.real(); 
  }

  if (m_coherentIntegralC) {
      complex_t nCC = sum_amps( m_normalisationsCC, m_C.matrixElements(), m_C.matrixElements() );
      nC = nCC.real(); 
  }

  mix = rAC * std::conj(sumFactor);




  real_t intTerm = 2 * mix.real();

  

  if (m_debug) INFO("interference = "<<intTerm);
  N = nA  + std::norm(sumFactor) * nC  + intTerm;

 // N = int_gp2().real() * nA + std::norm(sumFactor) * std::norm(qp) * int_gm2().real() * nC + 2 * std::real(int_gpm() * qp * mix);

  if (m_debug) INFO("Normalisation = "<<N);
  return N_slow;
}



/*
real_t tCoherentSum::norm() const {
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

  complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );


  real_t nA = m_A.norm();

  real_t nC = m_C.norm();


  if (m_coherentIntegralA) {
      complex_t nAA = sum_amps( m_normalisationsAA, m_A.matrixElements(), m_A.matrixElements() );
      nA = nAA.real(); 
  }
  
  if (m_coherentIntegralC) {
      complex_t nCC = sum_amps( m_normalisationsCC, m_C.matrixElements(), m_C.matrixElements() );
      nC = nCC.real(); 
  }
  
  complex_t mix = nAC *  std::conj(sumFactor);
  real_t intTerm = 2 * mix.real();

  

  if (m_debug) INFO("interference = "<<intTerm);
  real_t N = nA   + std::norm(sumFactor) * nC  + intTerm;

  if (m_debug) INFO("Normalisation = "<<N);
  return N;
}
*/


real_t tCoherentSum::testnorm() const {
  if (m_debug) INFO("Getting the value for the normalised pdf");
 
  if (m_debug) INFO("Getting SumFactor");
  complex_t sumFactor = getSumFactor();  
  if (m_debug) INFO("SumFactor = "<<sumFactor);
 if (m_debug) INFO("Get Matrix Elements for A,C");
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
    std::complex<double> rAC = 0; 
   auto eventsAC = m_integratorAC.events();

    complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );
  real_t nA = m_A.norm();
  if (m_debug) INFO("m_A norm = "<<nA);
  real_t nC = m_C.norm();
  if (m_debug) INFO("m_C norm = "<<nC);


  if (m_coherentIntegralA) {
      complex_t nAA = sum_amps( m_normalisationsAA, m_A.matrixElements(), m_A.matrixElements() );
      nA = nAA.real(); 
  }
  
  if (m_coherentIntegralC) {
      complex_t nCC = sum_amps( m_normalisationsCC, m_C.matrixElements(), m_C.matrixElements() );
      nC = nCC.real(); 
  }
  
//Note this part might crash - if so replace getVal with getValNoCache but I think its much slower!

    for( int i=0 ; i < eventsAC.size(); i++ ) {

      if (m_debug) INFO("Getting value for ab");
      if (m_debug) INFO("At "<<i<<" out of "<<eventsAC.size());
    if (m_debug) INFO("m_A matrixElements = "<<m_A.matrixElements().size()<<" long");
    if (m_debug) INFO("A = "<<m_A.getVal(eventsAC[i]));
    if (m_debug) INFO("C = "<<m_C.getVal(eventsAC[i]));

      auto ab = m_A.getValNoCache(eventsAC[i]) * std::conj(m_C.getValNoCache(eventsAC[i])) * exp(Constant(0,1)() * correction(eventsAC[i]))/(double)eventsAC.size();
  //    if (std::abs(ab) > 1e3){
    if (m_debug) INFO("ab = "<<ab);
 
    if (m_debug) INFO("Getting value for |a|");

    auto abA = abs(m_A.getValNoCache(eventsAC[i]));
    if (m_debug) INFO("|a| = "<<abA);
    if (m_debug) INFO("Getting value for |c|");
    auto abC = abs(m_C.getValNoCache(eventsAC[i]));

    if (m_debug) INFO("|c| = "<<abC);
    
      if (abA>1e10 || abC>1e10){
  //      INFO("A or C too large?");
      }
      else{
//   if (m_debug)   INFO("ab = "<<abs(ab));
//    if (m_debug)    INFO("correction = "<<correction(eventsAC[i]));
//      if (m_debug)    INFO("rAC = "<<abs(rAC));

        rAC = rAC + ab;
    }
    }


  complex_t mix = rAC  * std::conj(sumFactor);
  real_t intTerm = 2 * mix.real();

  

  if (m_debug) INFO("interference = "<<intTerm);
  real_t N = nA + std::norm(sumFactor) * nC  + intTerm;

  if (m_debug) INFO("Normalisation = "<<N);
  return N;
}

double tCoherentSum::slowNorm(){
  if (m_debug) INFO("Started slowNorm");
  if (m_debug) INFO("SumFactor prefix = "<<m_SFType);
  complex_t sumFactor = getSumFactor();  
  if (m_debug) INFO("sumFactor = "<<sumFactor);
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




   std::complex<double> rAC = 0; 
   if (m_debug) INFO("Getting events for AC");
   auto eventsAC = m_integratorAC.events();
   if (m_debug) INFO("Have "<<eventsAC.size()<< " integration events"); 
    for( int i=0 ; i < eventsAC.size(); i++ ) {
      if (m_debug) INFO("Getting value for ab");
      auto ab = m_A.getValNoCache(eventsAC[i]) * std::conj(m_C.getValNoCache(eventsAC[i])) * exp(Constant(0,1)() * correction(eventsAC[i]))/(double)eventsAC.size();
  //    if (std::abs(ab) > 1e3){
    if (m_debug) INFO("ab = "<<ab);
 
      if (m_debug) INFO("Getting value for |a|");

      auto abA = abs(m_A.getValNoCache(eventsAC[i]));
    if (m_debug) INFO("|a| = "<<abA);
      if (m_debug) INFO("Getting value for |c|");
      auto abC = abs(m_C.getValNoCache(eventsAC[i]));

    if (m_debug) INFO("|c| = "<<abC);
    
      if (abA>1e10 || abC>1e10){
  //      INFO("A or C too large?");
      }
      else{
//   if (m_debug)   INFO("ab = "<<abs(ab));
//    if (m_debug)    INFO("correction = "<<correction(eventsAC[i]));
//      if (m_debug)    INFO("rAC = "<<abs(rAC));

        rAC = rAC + ab;
    }
    }

    if (m_debug) INFO("AC = "<<rAC);



  if(m_debug) INFO("Do I have norm?");   
  double n_fast = m_norm;
  if (m_debug) INFO("m_norm = "<<n_fast);


//  double norm = m_A.norm() * m_B.norm() + m_C.norm() * m_D.norm();


  auto inter = rAC  * std::conj(sumFactor);
  auto Norm = m_Anorm  + m_Cnorm  ;


  Norm += 2 * std::real(inter);  
return Norm;
}

void tCoherentSum::reset(bool resetEvents){
  if (resetEvents){
    m_events1 = nullptr;

    m_integratorAC = Integrator<10>();

    m_integratorAA = Integrator<10>();

    m_integratorCC = Integrator<10>();

    m_A.reset(true);

    m_C.reset(true);

  }
}

void tCoherentSum::setEvents(EventList& list1){
  m_events1 = &list1;

  m_A.setEvents(list1);

  m_C.setEvents(list1);

}

void tCoherentSum::setMC(EventList& list1){
  m_A.setMC(list1);
  
  m_C.setMC(list1);
  m_integratorAC = Integrator<10>(&list1);
  
  m_integratorAA = Integrator<10>(&list1);
  
  m_integratorCC = Integrator<10>(&list1);
  

  
}

complex_t tCoherentSum::getVal(const Event& evt1) const {
  complex_t A = m_A.getValNoCache(evt1);

  complex_t C = m_C.getValNoCache(evt1);

  complex_t f = correction(evt1);
  
  auto i = Constant(0,1);
  auto corr = exp(i() * f);
  auto sumFactor = getSumFactor();
  complex_t val = A  *exp(i()*f/2.) + sumFactor* C  *exp(-i()*f/2.);
  auto param_x = m_mps["tCoherentSum::x"];
  auto param_y = m_mps["tCoherentSum::y"];



  auto param_absqp = m_mps["tCoherentSum::absqp"];

  auto param_phiqp = m_mps["tCoherentSum::phiqp"];




  auto D0Props = ParticleProperties::get(m_events1->eventType().mother());
  auto D0Width = D0Props->width();
  auto D0Lifetime = D0Props->lifetime();
  auto D0Mass = D0Props->mass();
  auto mix_x = param_x->mean();
  auto mix_y = param_y->mean();
  auto absqp = param_absqp->mean();
  auto phiqp = param_phiqp->mean();


  auto qp = absqp * fcn::exp(Constant(0, 1)() * phiqp)(); 
//  auto t = evt1[4 * m_events1->eventType().size()];
  auto t = evt1[4 * m_events1->eventType().size()];
  auto a = fcn::exp(-t/D0Lifetime - mix_y * t/(2*D0Lifetime) - Constant(0, 1) * D0Mass * t + Constant(0, 1) * mix_x * t/(2*D0Lifetime)); 
  auto b = fcn::exp(-1 * Constant(0,1) * mix_x * t/(2 * D0Lifetime) + mix_y * t/D0Lifetime);




//  auto a = fcn::exp(-D0Width * t - mix_y * D0Width * t/2 - i() * D0Mass * t + i * mix_x * D0Width * t/2); 
//  auto b = fcn::exp(-i() * mix_x * D0Width * t/2 + D0Width * mix_y * t);
 auto gplus = a * (1 + b);
auto  gminus = a * (1-b);

auto expr = A  *exp(i()*f/2.) * gplus + C *exp(-i()*f/2.)* gminus * sumFactor * qp;// + sumFactor* gminus * qp *C  *exp(-i()*f/2);
val = expr();
  //auto gplus = fcn::exp(-i() * )

  //INFO("Correction  = "<<f);
  if (m_debug) INFO("A2 = "<<std::norm(A));

  if (m_debug) INFO("C2 = "<<std::norm(C));

  if (m_debug) INFO("ReAC* = "<<std::real(A*std::conj(C)* corr * std::conj(sumFactor) ) );
  if (m_debug) INFO("val2 = "<<std::norm(val));
  if (m_debug) INFO("A2+C2-2ReAC* = "<<std::norm(A) + std::norm(C)*std::norm(sumFactor) - 2 *std::real(A*std::conj(C)*sumFactor*exp(i()*f)));

  if (m_flat){
    val = 1;
  }


  return val;
}
complex_t tCoherentSum::getValNoCache(const Event& evt1) const {
  complex_t A = m_A.getValNoCache(evt1);
  
  complex_t C = m_C.getValNoCache(evt1);
  
  complex_t f = correction(evt1);
  
  auto i = Constant(0,1);

  auto sumFactor = getSumFactor();
  complex_t val = A  *exp(i()*f/2.)  + sumFactor * C   *exp(-i()*f/2.);
  //complex_t val = A *exp(i()*f)* B - C * D;
  if (m_flat){
    val = 1;
  }
  return val;
}


void tCoherentSum::debugNorm()
{
auto i = Constant(0,1);
  auto n_slow = 0.0;
  auto nEvents = m_integratorAC.events().size();
  INFO("nEvents = "<<nEvents);
  for(unsigned int j = 0; j < m_integratorAC.events().size(); ++j ) {
    auto val = getValNoCache(m_integratorAC.events()[j]);
    auto val2 = val.real() * val.real() + val.imag() * val.imag();
    auto val3 = m_A.getValNoCache(m_integratorAC.events()[j]) * std::conj(m_C.getValNoCache(m_integratorAC.events()[j]));
    auto rval3 = m_A.getValNoCache(m_integratorAC.events()[j]) * std::conj(m_C.getValNoCache(m_integratorAC.events()[j]))/(std::abs(m_A.getValNoCache(m_integratorAC.events()[j])) * std::abs(m_C.getValNoCache(m_integratorAC.events()[j])));
    auto rpval3 = m_A.getValNoCache(m_integratorAC.events()[j]) * std::conj(m_C.getValNoCache(m_integratorAC.events()[j])) * exp(i() * correction(m_integratorAC.events()[j]))/(std::abs(m_A.getValNoCache(m_integratorAC.events()[j])) * std::abs(m_C.getValNoCache(m_integratorAC.events()[j])));
    auto pval3 = m_A.getValNoCache(m_integratorAC.events()[j]) * std::conj(m_C.getValNoCache(m_integratorAC.events()[j])) * exp(i() * correction(m_integratorAC.events()[j]));

    auto val5 = correction(m_integratorAC.events()[j]);
    auto cdd = std::real(rval3);
    auto dd = std::imag(log(val3/std::abs(val3)));
    auto cddp = std::real(rpval3);
    auto ddp = std::imag(log(pval3/std::abs(pval3)));
    auto diffdd = ddp - dd;
    n_slow += (double)val2/(double)nEvents;
    auto frac = 100*(double)j/(double)nEvents;
    if (m_debug){
      INFO("At "<<j <<" of "<<nEvents<<" events");
      INFO("ABC*D* = "<<val);
      INFO("|AB-CD|^2 = "<<val2);
      INFO("AC* = "<<val3);
      INFO("AC*e^(if) = "<<pval3);
      INFO("cos(dd) = "<<cdd); 
      INFO("dd = "<<dd); 
      INFO("cos(dd + f) = "<<cddp);
      INFO("dd + f = "<<ddp);
      INFO("dd + f - dd = "<<diffdd);
      INFO("n_slow = "<<n_slow);

      INFO("exp(i f(s+,s-)="<<val5);
    }
  }
  auto integrate = []( const EventList& events, const CoherentSum& A, const CoherentSum& B, bool debug=false) -> double
  {
   auto r = 0.0; 
   auto N = events.size(); 
    for( int i=0 ; i < N; i++ ) {
     auto ab = A.getValNoCache(events[i]) * std::conj(B.getValNoCache(events[i]));


     auto rab = std::real(ab)/N;


  long double L=events.size();
    
    
     r = r + rab;

   if (std::abs(rab)<L )  r = r + rab;

    }

    return r;
  };
  double nAA = std::real( integrate( m_integratorAC.events(), m_A, m_A, m_debug ) );

  double nCC = std::real( integrate( m_integratorAC.events(), m_C, m_C, m_debug ) );

  auto   nAC = ( integrate( m_integratorAC.events(), m_A, m_C, m_debug ) );

  
  INFO("integral |A|^2 = "<<nAA);
  INFO("fast |A|^2 = "<<m_A.norm());


  INFO("|C| = "<<nCC);
  INFO("fast |C|^2 = "<<m_C.norm());


  INFO("AC* = "<<-2*nAC);
  double I = ( nAA  + nCC  - 2*std::real(nAC) );
  double n_fast = m_norm;
  auto n_slow2 = slowNorm();
  double difference = std::abs(n_fast - I)/n_fast * 100.;
  INFO("My value for n_fast = "<<n_fast);
  INFO("My value for I = "<<I);
  INFO("My value for n_slow = "<<n_slow2);
  INFO("The difference = "<<difference<<"%");
}


void tCoherentSum::debug(const Event& evt1) const 
{
  INFO( "A[x, y] = " << getValNoCache(evt1) << " " << getValNoCache(evt1) << " " << m_norm << " " << prob(evt1) << " " << prob_unnormalised(evt1) );
  INFO("|AB-CD|^2 = "<<prob(evt1));
  INFO("|A|^2 = "<<m_A.prob(evt1));

  INFO("|C|^2 = "<<m_C.prob(evt1));

  INFO("Norm = "<<m_norm);

  if (m_debug) INFO("Getting the value for the normalised pdf");
   auto eventsAC = m_integratorAC.events();
  if (m_debug) INFO("Get sumFactor");
 
  complex_t sumFactor = getSumFactor();  
  INFO("sumFactor="<<sumFactor);

  INFO("Get TD params");
  INFO("Get x");
   auto param_x = m_mps["tCoherentSum::x"];

  INFO("x = "<<param_x->mean());
  auto param_y = m_mps["tCoherentSum::y"];
  INFO("y = "<<param_y->mean());
  auto param_absqp = m_mps["tCoherentSum::absqp"];
  INFO("|q/p| = "<<param_absqp->mean());
  auto param_phiqp = m_mps["tCoherentSum::phiqp"];
  INFO("phi(q/p) = "<<param_phiqp->mean());




  auto evtType = eventsAC.eventType();
  INFO("eventType "<<evtType);
  INFO("Get mother "<<eventsAC.eventType().mother()<<" properties");

  auto D0Props = ParticleProperties::get(eventsAC.eventType().mother());
  auto D0Width = D0Props->width();
  auto D0Lifetime = D0Props->lifetime();
  INFO(eventsAC.eventType().mother()<<" width = "<<D0Width);
  
  auto D0Mass = D0Props->mass();

auto h = 6.62607015e-34;
auto c = 299792458.0;
D0Mass = D0Mass * pow(c, 2)/h;
  INFO(eventsAC.eventType().mother()<<" mass = "<<D0Mass);
  auto mix_x = param_x->mean();
  auto mix_y = param_y->mean();
  auto absqp = param_absqp->mean();
  auto phiqp = param_phiqp->mean();
  auto qp = absqp * fcn::exp( Constant(0, 1)() * phiqp)(); 
//  auto a = fcn::exp(-D0Width * t - mix_y * D0Width * t/2 - i() * D0Mass * t + i * mix_x * D0Width * t/2); 
//  auto b = fcn::exp(-i() * mix_x * D0Width * t/2 + D0Width * mix_y * t);
//  auto gplus = a * (1 + b);
// auto  gminus = a * (1-b);  


  INFO("Calculate TD parts of normalisation");
//  auto int_gp2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  + 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gp2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) + 1/(1 + fcn::pow(mix_x, 2)));


  //auto int_gm2 = 1/(2/D0Lifetime + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  - 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );

  //auto int_gm2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * (D0Lifetime/(1 - fcn::pow(mix_y, 2)) - D0Lifetime/(1 + fcn::pow(mix_x, 2))); //  (2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) - 1/(1 + fcn::pow(mix_x, 2)));

  //auto int_gpm = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime) - Constant(0,1)() * mix_x * (1/D0Lifetime)/(4 * fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x/D0Lifetime, 2)()/4);

  //auto int_gp2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) + 1/(1 + fcn::pow(mix_x, 2)));
  auto int_gp2 = (1/2)  * (1/(1 - fcn::pow(mix_y, 2)) + 1/(1 + fcn::pow(mix_x, 2)));
  INFO("Int |g_+(t)|^2 dt = "<<int_gp2());

  //auto int_gm2 = 1/(2/D0Lifetime + mix_y * D0Width) + 1/(-2 * D0Width + mix_y * D0Width)  - 4 * D0Width/(4*fcn::pow(D0Width, 2)() + fcn::pow(mix_x * D0Width, 2)/4 );

  //auto int_gm2 = 1/(2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * (D0Lifetime/(1 - fcn::pow(mix_y, 2)) - D0Lifetime/(1 + fcn::pow(mix_x, 2))); //  (2/D0Lifetime + mix_y/D0Lifetime) + 1/(-2/D0Lifetime + mix_y/D0Lifetime)  - 4 * (1/D0Lifetime)/(4*fcn::pow(1/D0Lifetime, 2)() + fcn::pow(mix_x / D0Lifetime, 2)/4 );
  //auto int_gm2 = (1/2) * D0Lifetime * (1/(1 - fcn::pow(mix_y, 2)) - 1/(1 + fcn::pow(mix_x, 2)));
  auto int_gm2 = (1/2) * (1/(1 - fcn::pow(mix_y, 2)) - 1/(1 + fcn::pow(mix_x, 2)));
  INFO("Int |g_-(t)|^2 dt = "<<int_gp2());
  //auto int_gpm = D0Lifetime/(2 + mix_y) + D0Lifetime/(-2 + mix_y) - Constant(0,1)() * mix_x * D0Lifetime/(4 + fcn::pow(mix_x, 2)()/4);
  auto int_gpm = 1/(2 + mix_y) + 1/(-2 + mix_y) - Constant(0,1)() * mix_x * 1/(4 + fcn::pow(mix_x, 2)()/4);
  INFO("Int g_+(t)g_-^*(t) dt = "<<int_gpm());



 INFO("Get norms for A,C");



  real_t nA = m_A.norm();

  real_t nC = m_C.norm();


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
  INFO("Get interference without corrections");
  complex_t nAC = sum_amps( m_normalisationsAC, m_A.matrixElements(), m_C.matrixElements() );

  INFO("interference = "<<nAC);

  complex_t mix = nAC * std::conj(sumFactor);
  if (m_debug) INFO("Calculate uncorrected normalisation");
  real_t N = int_gp2().real() * nA + std::norm(sumFactor) * std::norm(qp) * int_gm2().real() * nC + 2 * std::real(int_gpm() * qp * mix);

  INFO("Norm = "<<N);

  real_t NAC = nA + nC + 2*nAC.real();
  INFO("nA = "<<nA);
  INFO("nC = "<<nC);
  INFO("N(A+C) = "<<NAC);

  auto t0 = evt1[4 * m_events1->eventType().size()];
  auto a0 = fcn::exp(- t0/D0Lifetime - mix_y * t0/(2*D0Lifetime) - Constant(0, 1) * D0Mass * t0 + Constant(0, 1) * mix_x * t0/(2*D0Lifetime)); 
  auto b0 = fcn::exp(-1 * Constant(0,1) * mix_x * t0/(2 * D0Lifetime) + mix_y * t0/D0Lifetime);
 auto gplus0 = a0 * (1 + b0);
auto  gminus0 = a0 * (1-b0);

auto slow_gp2_int = 0;
auto slow_gm2_int = 0;
complex_t slow_gpm_int = 0;
double N_slow =0;
for (int i=0 ; i<eventsAC.size() ; i++){
  auto t = eventsAC[i][4 * m_events1->eventType().size()];
  auto a = fcn::exp(- t/D0Lifetime - mix_y * t/(2*D0Lifetime) - Constant(0, 1) * D0Mass * t + Constant(0, 1) * mix_x * t/(2*D0Lifetime)); 
  auto b = fcn::exp(-1 * Constant(0,1) * mix_x * t/(D0Lifetime) + mix_y * t/D0Lifetime);
 auto gplus = a * (1 + b);
auto  gminus = a * (1-b);

  slow_gp2_int += std::norm(gplus());
  slow_gm2_int += std::norm(gminus());
  slow_gpm_int += gplus() * std::conj(gminus());
    auto v = std::norm(getVal(eventsAC[i]));
    N_slow += v;

}



INFO("time = "<<t0);
INFO("a0 = "<<a0());
INFO("b0 = "<<b0());
INFO("g_+(t)0 = "<<gplus0());
INFO("g_-(t)0 = "<<gminus0());

INFO("Int |g_+|^2 dt = "<<slow_gp2_int);
INFO("Int |g_-|^2 dt = "<<slow_gp2_int);
INFO("Int g_+g_-^* dt = "<<slow_gpm_int);
INFO("nEvents = "<<eventsAC.size());
double nevents = eventsAC.size();
INFO("Int |g_+|^2 dt = "<<slow_gp2_int/nevents);
INFO("Int |g_-|^2 dt = "<<slow_gp2_int/nevents);
INFO("Int g_+g_-^* dt = "<<slow_gpm_int/nevents);


auto fast_gp2_int = (1/2.) * D0Lifetime * (1/(1-pow(mix_y, 2)) + 1/(1 - pow(mix_x, 2)));
auto fast_gm2_int = std::norm(qp) * (1/2.) *  D0Lifetime * (1/(1-pow(mix_y, 2)) - 1/(1 - pow(mix_x, 2)));
INFO("D0Lifetime = "<<D0Lifetime);
INFO("Fast |g_+|^2 dt = "<<fast_gp2_int);
INFO("Fast |g_-|^2 dt = "<<fast_gm2_int);
INFO("N_slow = "<<N_slow);
}

std::vector<std::vector<FitFraction> > tCoherentSum::fitFractions(const LinearErrorPropagator& linProp){
    outputFractionsA = m_A.fitFractions(linProp);
    
    outputFractionsC = m_C.fitFractions(linProp);
    
    std::vector<std::vector<FitFraction> > outputFractions;   
    outputFractions.push_back(outputFractionsA);
    
    outputFractions.push_back(outputFractionsC);
    
    return outputFractions;
}
