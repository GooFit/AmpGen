#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
using namespace std;
using namespace AmpGen;

CorrelatedSum::CorrelatedSum(){
    
}
CorrelatedSum::CorrelatedSum(const EventType& type1, const EventType& type2, const MinuitParameterSet& mps):
    m_A(type1, mps),
    m_B(type2.conj(true), mps),
    m_C(type1.conj(true), mps),
    m_D(type2, mps),
    m_debug(NamedParameter<bool>("CorrelatedSum::debug", false, "Print Debug messages for CorrelatedSum")),
    m_coherentIntegral(NamedParameter<bool>("CorrelatedSum::coherentIntegral", false, "Do the integral for each coherent sum"))
    {
        
       
         m_normalisationsAC.resize( m_A.matrixElements().size(), m_C.matrixElements().size() ); 
         m_normalisationsBD.resize( m_B.matrixElements().size(), m_D.matrixElements().size() ); 

         m_normalisationsAA.resize( m_A.matrixElements().size(), m_A.matrixElements().size() ); 
         m_normalisationsBB.resize( m_B.matrixElements().size(), m_B.matrixElements().size() ); 
         m_normalisationsCC.resize( m_C.matrixElements().size(), m_C.matrixElements().size() ); 
         m_normalisationsDD.resize( m_D.matrixElements().size(), m_D.matrixElements().size() ); 

    }
void CorrelatedSum::prepare(){
        m_A.prepare();
        m_B.prepare();
        m_C.prepare();
        m_D.prepare();
 
    if (m_debug) INFO("Preparing CorrelatedSum");
    std::vector<unsigned int> changedPdfIndicesA;
    std::vector<unsigned int> changedPdfIndicesB;
    std::vector<unsigned int> changedPdfIndicesC;
    std::vector<unsigned int> changedPdfIndicesD;
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
    for ( unsigned int i = 0; i < matrixElementsA.size(); ++i ) {
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
        //else if ( i == 0 && m_verbosity ) WARNING( "No data events specified for " << this );
        clockThisElement.stop();
        /*
        if ( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ) {
        if (m_debug) INFO( amp.name() << " (t = " << clockThisElement << " ms, nCalls = " << m_prepareCalls << ", events = " << m_events->size() << ")" );
        printed = true;
        }
        */
        changedPdfIndicesA.push_back(i);
        amp.resetExternals();
        }

    if (m_debug) INFO("Get changedPdfIndiciesB for "<<matrixElementsB.size()<<" Elements");

    for ( unsigned int i = 0; i < matrixElementsB.size(); ++i ) {
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
        //else if ( i == 0 && m_verbosity ) WARNING( "No data events specified for " << this );
        clockThisElement.stop();
        /*
        if ( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ) {
        if (m_debug) INFO( amp.name() << " (t = " << clockThisElement << " ms, nCalls = " << m_prepareCalls << ", events = " << m_events->size() << ")" );
        printed = true;
        }
        */
        changedPdfIndicesB.push_back(i);
        amp.resetExternals();
        }

    if (m_debug) INFO("Get changedPdfIndiciesC");
    for ( unsigned int i = 0; i < matrixElementsC.size(); ++i ) {
        auto& amp = matrixElementsC[i].amp;
        amp.prepare();
        if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
        ProfileClock clockThisElement;
        if ( m_events1 != nullptr ) {
            if ( matrixElementsC[i].addressData == 999 ) 
                matrixElementsC[i].addressData = m_events1->registerExpression( amp );
            m_events1->updateCache( amp, matrixElementsC[i].addressData );
            } 
        //else if ( i == 0 && m_verbosity ) WARNING( "No data events specified for " << this );
        clockThisElement.stop();
        /*
        if ( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ) {
        if (m_debug) INFO( amp.name() << " (t = " << clockThisElement << " ms, nCalls = " << m_prepareCalls << ", events = " << m_events->size() << ")" );
        printed = true;
        }
        */
        changedPdfIndicesC.push_back(i);
        amp.resetExternals();
        }
        
    if (m_debug) INFO("Get changedPdfIndiciesD");
    for ( unsigned int i = 0; i < matrixElementsD.size(); ++i ) {
        auto& amp = matrixElementsD[i].amp;
        amp.prepare();
        if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
        ProfileClock clockThisElement;
        if ( m_events2 != nullptr ) {
        if ( matrixElementsD[i].addressData == 999 ) 
            matrixElementsD[i].addressData = m_events2->registerExpression( amp );
        m_events2->updateCache( amp, matrixElementsD[i].addressData );
        } 
        //else if ( i == 0 && m_verbosity ) WARNING( "No data events specified for " << this );
        clockThisElement.stop();
        /*
        if ( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ) {
        if (m_debug) INFO( amp.name() << " (t = " << clockThisElement << " ms, nCalls = " << m_prepareCalls << ", events = " << m_events->size() << ")" );
        printed = true;
        }
        */
        changedPdfIndicesD.push_back(i);
        amp.resetExternals();
        }
        if (m_debug) INFO("Get Normalisation");
        //if (m_integratorAC.isReady() && m_integratorBC.isReady()){
        updateNorms(changedPdfIndicesA, changedPdfIndicesB, changedPdfIndicesC, changedPdfIndicesD);
        //}
        m_norm = norm();
        
        m_prepareCalls++;
}

void CorrelatedSum::updateNorms(const std::vector<unsigned int>& iA, const std::vector<unsigned int>& iB,
                                const std::vector<unsigned int>& iC, const std::vector<unsigned int>& iD){
    
    //Updating Norms - adapted from CoherentSum where we need to perform different integrals to the CS
    /*
    for (auto& i : iA ) m_integratorAA.prepareExpression(m_A.matrixElements()[i].pdf);
    for (auto& i : iB ) m_integratorBB.prepareExpression(m_B.matrixElements()[i].pdf);
    for (auto& i : iC ) m_integratorCC.prepareExpression(m_C.matrixElements()[i].pdf);
    for (auto& i : iD ) m_integratorDD.prepareExpression(m_D.matrixElements()[i].pdf);
    */


    for (auto& i : iA ) m_integratorAC.prepareExpression(m_A.matrixElements()[i].amp);
    for (auto& i : iB ) m_integratorBD.prepareExpression(m_B.matrixElements()[i].amp);
    for (auto& i : iC ) m_integratorAC.prepareExpression(m_C.matrixElements()[i].amp);
    for (auto& i : iD ) m_integratorBD.prepareExpression(m_D.matrixElements()[i].amp);

    for (auto& i : iA ) m_integratorAA.prepareExpression(m_A.matrixElements()[i].amp);
    for (auto& i : iB ) m_integratorBB.prepareExpression(m_B.matrixElements()[i].amp);
    for (auto& i : iC ) m_integratorCC.prepareExpression(m_C.matrixElements()[i].amp);
    for (auto& i : iD ) m_integratorDD.prepareExpression(m_D.matrixElements()[i].amp);

    std::vector<size_t> cacheIndexA;
    std::vector<size_t> cacheIndexB;
    std::vector<size_t> cacheIndexC;
    std::vector<size_t> cacheIndexD;

    if (m_debug) INFO("integratorAC has "<<m_integratorAC.events().size());
    if (m_debug) INFO("integratorBD has "<<m_integratorBD.events().size());

/*
    if (m_debug) INFO("Print the pdf names that getCacheIndex is being fed");
    for (auto element : m_A.matrixElements()){
        if (m_debug) INFO("pdf name = "<<element.pdf.name());
        if (m_debug) INFO("getCacheIndex output = "<<m_integratorAC.events().getCacheIndex(element.pdf));
        cacheIndexA.push_back(m_integratorAC.events().getCacheIndex(element.pdf));
    }
    if (m_debug) INFO("Print the pdf names that getCacheIndex is being fed");
    for (auto element : m_B.matrixElements()){
        if (m_debug) INFO("pdf name = "<<element.pdf.name());
        if (m_debug) INFO("getCacheIndex output = "<<m_integratorBD.events().getCacheIndex(element.pdf));
        cacheIndexB.push_back(m_integratorBD.events().getCacheIndex(element.pdf));
    }
    if (m_debug) INFO("Print the pdf names that getCacheIndex is being fed");
    for (auto element : m_C.matrixElements()){
        if (m_debug) INFO("pdf name = "<<element.pdf.name());
        if (m_debug) INFO("getCacheIndex output = "<<m_integratorAC.events().getCacheIndex(element.pdf));
        cacheIndexC.push_back(m_integratorAC.events().getCacheIndex(element.pdf));
    }
    if (m_debug) INFO("Print the pdf names that getCacheIndex is being fed");
    for (auto element : m_D.matrixElements()){
        if (m_debug) INFO("pdf name = "<<element.pdf.name());
        if (m_debug) INFO("getCacheIndex output = "<<m_integratorBD.events().getCacheIndex(element.pdf));
        cacheIndexD.push_back(m_integratorBD.events().getCacheIndex(element.pdf));
    }
    */
    
    if (m_debug) INFO("Beginning to get CacheIndex for amplitudes A,B,C,D - for the integrator to use");
    if (m_debug) INFO("A index");
   /* 
    std::transform( m_A.matrixElements().begin(), m_A.matrixElements().end(), std::back_inserter(cacheIndexA) ,
        [this](auto& m) {return this->m_integratorAC.events().getCacheIndex( m.pdf ); } );
    
    if (m_debug) INFO("B index");
    std::transform( m_B.matrixElements().begin(), m_B.matrixElements().end(), std::back_inserter(cacheIndexB) ,
        [this](auto& m) {return this->m_integratorBD.events().getCacheIndex( m.pdf ); } );

    if (m_debug) INFO("C index");
    std::transform( m_C.matrixElements().begin(), m_C.matrixElements().end(), std::back_inserter(cacheIndexC) ,
        [this](auto& m) {return this->m_integratorAC.events().getCacheIndex( m.pdf ); } );

    if (m_debug) INFO("D index");
    std::transform( m_D.matrixElements().begin(), m_D.matrixElements().end(), std::back_inserter(cacheIndexD) ,
        [this](auto& m) {return this->m_integratorBD.events().getCacheIndex( m.pdf ); } );
    
    if (m_debug) INFO("Queuing the integral for nAC"); 
    for (auto i : iA ){
        for (auto j : iC){
            if (m_debug) INFO("At ("<<i<<", "<<j<<") for calculation using m_normalisationsAC" << &m_normalisationsAC);
            m_integratorAC.queueIntegral(cacheIndexA[i], cacheIndexC[j], i, j, &m_normalisationsAC);
        }
    }
    */


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
            m_integratorAC.queueIntegral(cacheIndexA[i], cacheIndexC[j], i, j, &m_normalisationsAC);
        }
    }
    if (m_debug) INFO("Finishing integral for nAC");
    m_integratorAC.flush();
    m_normalisationsAC.resetCalculateFlags();



   if (m_debug) INFO("Queuing the integral for nBD"); 
   for (auto i : iB ){
        for (auto j : iD){
            m_integratorBD.queueIntegral(cacheIndexB[i], cacheIndexD[j], i, j, &m_normalisationsBD);
        }
    }

    if (m_debug) INFO("Finishing integral for nBD");
    m_integratorBD.flush();
    m_normalisationsBD.resetCalculateFlags();

    if (m_coherentIntegral){

        if (m_debug) INFO("Queuing the integral for nAA"); 
        for (auto i : iA ){
            for (auto j : iA){
                m_integratorAA.queueIntegral(cacheIndexA[i], cacheIndexA[j], i, j, &m_normalisationsAA);
            }
        }
        if (m_debug) INFO("Finishing integral for nAA");
        m_integratorAA.flush();
        m_normalisationsAA.resetCalculateFlags();

        if (m_debug) INFO("Queuing the integral for nBB"); 
        for (auto i : iB ){
            for (auto j : iB){
                m_integratorBB.queueIntegral(cacheIndexB[i], cacheIndexB[j], i, j, &m_normalisationsBB);
            }
        }
        if (m_debug) INFO("Finishing integral for nBB");
        m_integratorBB.flush();
        m_normalisationsBB.resetCalculateFlags();

        if (m_debug) INFO("Queuing the integral for nCC"); 
        for (auto i : iC ){
            for (auto j : iC){
                m_integratorCC.queueIntegral(cacheIndexC[i], cacheIndexC[j], i, j, &m_normalisationsCC);
            }
        }
        if (m_debug) INFO("Finishing integral for nCC");
        m_integratorCC.flush();
        m_normalisationsCC.resetCalculateFlags();

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

real_t CorrelatedSum::norm() const {
    if (m_debug) INFO("Getting the value for the normalised pdf");
    complex_t nAC(0,0);
    complex_t nBD(0,0);

    complex_t nAA(0,0);
    complex_t nBB(0,0);
    complex_t nCC(0,0);
    complex_t nDD(0,0);
    //complex_t intTerm(0,0);
    if (m_debug) INFO("Get Matrix Elements for A,B,C,D");

    auto matrixA = m_A.matrixElements();
    auto matrixB = m_B.matrixElements();
    auto matrixC = m_C.matrixElements();
    auto matrixD = m_D.matrixElements();
    for (size_t i=0; i<matrixA.size(); i++){
        for (size_t j=0; j<matrixC.size(); j++){
            auto val = m_normalisationsAC.get(i, j) * matrixA[i].coefficient * std::conj(matrixC[j].coefficient);
            nAC += val;
        }
    }
    for (size_t i=0; i<matrixB.size(); i++){
        for (size_t j=0; j<matrixD.size(); j++){
            auto val = m_normalisationsBD.get(i, j) * matrixB[i].coefficient * std::conj(matrixD[j].coefficient);
            nBD += val;
        }
    }

    if (m_coherentIntegral){
        for (size_t i=0; i<matrixA.size(); i++){
            for (size_t j=0; j<matrixA.size(); j++){
                auto val = m_normalisationsAA.get(i, j) * matrixA[i].coefficient * std::conj(matrixA[j].coefficient);
                nAA += val;
            }
        }
        for (size_t i=0; i<matrixB.size(); i++){
            for (size_t j=0; j<matrixB.size(); j++){
                auto val = m_normalisationsBB.get(i, j) * matrixB[i].coefficient * std::conj(matrixB[j].coefficient);
                nBB += val;
            }
        }
        for (size_t i=0; i<matrixC.size(); i++){
            for (size_t j=0; j<matrixC.size(); j++){
                auto val = m_normalisationsCC.get(i, j) * matrixC[i].coefficient * std::conj(matrixC[j].coefficient);
                nCC += val;
            }
        }
        for (size_t i=0; i<matrixD.size(); i++){
            for (size_t j=0; j<matrixD.size(); j++){
                auto val = m_normalisationsDD.get(i, j) * matrixD[i].coefficient * std::conj(matrixD[j].coefficient);
                nDD += val;
            }
        }
    }
    else{
        nAA = m_A.norm();
        nBB = m_B.norm();
        nCC = m_C.norm();
        nDD = m_D.norm();
    }

    real_t nA = nAA.real(); 
    real_t nB = nBB.real();
    real_t nC = nCC.real();
    real_t nD = nDD.real();

    if (m_debug) INFO("Each normalisation is" );
    if (m_debug) INFO("A = "<<nAA);
    if (m_debug) INFO("B = "<<nBB);
    if (m_debug) INFO("C = "<<nCC);
    if (m_debug) INFO("D = "<<nDD);
    if (m_debug) INFO("AC = "<<nAC);
    if (m_debug) INFO("BD = "<<nBD);
    
    if (m_debug) INFO("Comparing individual normalisations");
    if (m_debug) INFO("calcA = "<<nA<<" m_A.norm() = "<<m_A.norm());
    if (m_debug) INFO("calcB = "<<nB<<" m_B.norm() = "<<m_B.norm());
    if (m_debug) INFO("calcC = "<<nC<<" m_C.norm() = "<<m_C.norm());
    if (m_debug) INFO("calcD = "<<nD<<" m_D.norm() = "<<m_D.norm());
   

    real_t intTerm = -2 * (nAC * nBD).real();
    if (m_debug) INFO("interference = "<<intTerm);
    real_t N = nA * nB + nC * nD + intTerm;
    if (m_debug) INFO("Normalisation = "<<N);
    return N;
}

void CorrelatedSum::reset(bool resetEvents){
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

void CorrelatedSum::setEvents(EventList& list1, EventList& list2){
    m_events1 = &list1;
    m_events2 = &list2;
    m_A.setEvents(list1);
    m_B.setEvents(list2);
    m_C.setEvents(list1);
    m_D.setEvents(list2);
}


void CorrelatedSum::setMC(EventList& list1, EventList& list2){
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

complex_t CorrelatedSum::getVal(const Event& evt1, const Event& evt2) const {
    complex_t A = m_A.getVal(evt1);
    complex_t B = m_B.getVal(evt2);
    complex_t C = m_C.getVal(evt1);
    complex_t D = m_D.getVal(evt2);
    complex_t val = A * B - C * D;
    return val;
}
complex_t CorrelatedSum::getValNoCache(const Event& evt1, const Event& evt2) const {
    complex_t A = m_A.getValNoCache(evt1);
    complex_t B = m_B.getValNoCache(evt2);
    complex_t C = m_C.getValNoCache(evt1);
    complex_t D = m_D.getValNoCache(evt2);
    complex_t val = A * B - C * D;
    return val;
}


void CorrelatedSum::debugNorm(){
     INFO("Testing the interference parameter in the correlated sum - we assume that the coherent sum's are actually normalised properly!");
    complex_t nAC_slow(0,0);
    complex_t nBD_slow(0,0);
    complex_t nAA_slow(0,0);
    complex_t nBB_slow(0,0);
    complex_t nCC_slow(0,0);
    complex_t nDD_slow(0,0);


    INFO("Initial value for nAC = "<<nAC_slow);
    INFO("Initial value for nBD = "<<nBD_slow);
    INFO("Initial value for nAA = "<<nAA_slow);
    INFO("Initial value for nBB = "<<nBB_slow);
    INFO("Initial value for nCC = "<<nCC_slow);
    INFO("Initial value for nDD = "<<nDD_slow);

    double i=0.0;

    for (auto& evt : m_integratorAC.events()){
        // if (m_debug) INFO("A = "<<m_A.getVal(evt));
        // if (m_debug) INFO("C = "<<m_C.getVal(evt));
       //  if (m_debug) INFO("slow = "<<evt.weight() * m_A.getVal(evt) * std::conj(m_C.getVal(evt))/evt.genPdf());
        nAC_slow += evt.weight() * m_A.getVal(evt) * std::conj(m_C.getVal(evt))/evt.genPdf();
        nAA_slow += evt.weight() * m_A.getVal(evt) * std::conj(m_A.getVal(evt))/evt.genPdf();
        nCC_slow += evt.weight() * m_C.getVal(evt) * std::conj(m_C.getVal(evt))/evt.genPdf();

        if (m_debug) INFO("nAA_slow = "<<nAA_slow);
        if (m_debug) INFO("nCC_slow = "<<nCC_slow);
        if (m_debug) INFO("nAC_slow = "<<nAC_slow);
        i++;
    }
    nAC_slow /= i;
    nAA_slow /= i;
    nCC_slow /= i;

    i=0;
    for (auto& evt : m_integratorBD.events()){
/*         if (m_debug) INFO("B = "<<m_B.getVal(evt));
         if (m_debug) INFO("Bconj = "<<std::conj(m_B.getVal(evt)));
         if (m_debug) INFO("D = "<<m_D.getVal(evt));
         if (m_debug) INFO("Dconj = "<<std::conj(m_D.getVal(evt)));
        if (m_debug) INFO("nBB_slow = "<<m_B.getVal(evt) * std::conj(m_B.getVal(evt)) << " \t sum  =" <<nBB_slow);
        if (m_debug) INFO("nDD_slow = "<<m_D.getVal(evt) * std::conj(m_D.getVal(evt)) << " \t sum  =" <<nDD_slow);
        if (m_debug) INFO("nBD_slow = "<<m_B.getVal(evt) * std::conj(m_D.getVal(evt)) << " \t sum  =" <<nBD_slow);
 */
        nBD_slow += evt.weight() * m_B.getVal(evt) * std::conj(m_D.getVal(evt))/evt.genPdf();
        nBB_slow += evt.weight() * m_B.getVal(evt) * std::conj(m_B.getVal(evt))/evt.genPdf();
        nDD_slow += evt.weight() * m_D.getVal(evt) * std::conj(m_D.getVal(evt))/evt.genPdf();
        if (m_debug) INFO("nBB_slow = "<<nBB_slow);
        if (m_debug) INFO("nDD_slow = "<<nDD_slow);
        if (m_debug) INFO("nBD_slow = "<<nBD_slow);

        
        
        i++;
    }
    INFO("i = "<<i);
    nBD_slow /= i;
    nBB_slow /= i;
    nDD_slow /= i;
    {
         INFO("slow A = "<<nAA_slow.real());
         INFO("fast A = "<<m_A.norm());
         INFO("slow B = "<<nBB_slow.real());
         INFO("fast B = "<<m_B.norm());
         INFO("slow C = "<<nCC_slow.real());
         INFO("fast C = "<<m_C.norm());
         INFO("slow D = "<<nDD_slow.real());
         INFO("fast D = "<<m_D.norm());
         INFO("nBD_slow = "<<nBD_slow);
    }
    //auto nACBD_slow = nAC_slow * nBD_slow;
    real_t inter_slow = -2 * (nAC_slow * nBD_slow).real();
    real_t inter_fast = m_norm - m_A.norm() * m_B.norm() - m_C.norm() * m_D.norm();
    {
         INFO("nACBD_slow = "<<inter_slow);
         INFO("nACBD_fast = "<<inter_fast);
    }

    real_t n_slow = nAA_slow.real() * nBB_slow.real() + nCC_slow.real() * nDD_slow.real() + inter_slow;
    double n_fast = m_norm;
    double difference = std::abs(n_fast - n_slow)/n_fast * 100.;
    {
         INFO("My value for n_slow = "<<n_slow);
         INFO("My value for n_fast = "<<n_fast);
         INFO("The difference = "<<difference<<"%");
    }
}


