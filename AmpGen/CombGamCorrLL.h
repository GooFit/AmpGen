#ifndef COMBGAMCORRLL
#define COMBGAMCORRLL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/CorrelatedLL.h"
#include "TRandom3.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PhaseCorrection.h"

#include <tuple>
namespace AmpGen
{
  class EventList;
  /**
   * @class CombLL
   * @brief A combined log-likelihood for elated amplitudes.
   **/
  class CombGamCorrLL
  {
      private:
        EventType m_SigType;
        std::vector<EventType> m_TagType;
        MinuitParameterSet m_mps;
        bool m_debug;
        std::vector<std::vector<complex_t> >m_A={};
        std::vector<std::vector<complex_t> >m_B={};
        std::vector<std::vector<complex_t> >m_C={};
        std::vector<std::vector<complex_t> >m_D={};
        std::vector<std::vector<complex_t> >m_gA={};
        std::vector<std::vector<complex_t> >m_gC={};

        real_t m_ANorm=0;
        std::map<size_t, real_t> m_BNorms={};
        real_t m_CNorm=0;
        std::map<size_t, real_t> m_DNorms={};

        std::vector<complex_t> m_AMC={};
        std::map<size_t, std::vector<complex_t> > m_BMC={};
        std::vector<complex_t> m_CMC={};
        std::map<size_t, std::vector<complex_t> >m_DMC={};

        std::vector<bool> m_constAmp;


        std::vector<int> m_gammaSigns;
        std::vector<int> m_useXYs;
        std::vector<int> m_BConj;
        bool m_print;

        std::vector<PhaseCorrection> m_pc={};
        std::vector<PhaseCorrection> m_pcg = {};
        PhaseCorrection m_pcMC;
        PhaseCorrection m_pcT;
        


      public:
        CombGamCorrLL() = default;
        CombGamCorrLL(std::vector<EventList> SigData, 
		   std::vector<EventList> TagData, 

		   EventType SigType,
		   std::vector<EventType> TagType,
           std::vector<bool> constAmp,

		   std::vector<EventList> GamData, 
           std::vector<int> gammaSigns,
           std::vector<int> useXYs,

           std::vector<int> BConj,
            MinuitParameterSet mps):
            m_SigType(SigType),
            m_TagType(TagType),
        
            m_mps(mps),
            m_constAmp(constAmp),

            m_gammaSigns(gammaSigns),
            m_useXYs(useXYs),
            m_BConj(BConj),
            m_debug(NamedParameter<bool>("CombGamCorrLL::Debug", false, "Debug CombLL")),
            m_print(NamedParameter<bool>("CombGamCorrLL::Print", false, "Print CombLL"))
			
            {
                //MC for AC
            TRandom3 rndm;
            rndm.SetSeed( NamedParameter<size_t>("CombGamCorrLL::Seed", 0) );
            gRandom = &rndm;
            size_t NInt = NamedParameter<size_t>("CombGamCorrLL::NInt", 1e5);
            EventList mc =  Generator<>(SigType, &rndm).generate(NInt);
            INFO("Make phase correction for MC using "<<NInt<<" integration events");
            m_pcMC = PhaseCorrection(m_mps);
            m_pcMC.setEvents(mc);
            m_pcMC.prepareCache();
            INFO("Prepared phase correction");
            INFO("Setting up MC for A");
            CoherentSum * A = new CoherentSum(m_SigType, m_mps);            A->setEvents(mc); A->setMC(mc); A->transferParameters(); A->prepare();  m_ANorm = A->norm(); for (auto& evt:mc) {m_AMC.push_back((complex_t)A->getValNoCache(evt));} A->reset(true); delete A;
           // for (auto v : m_AMC) m_ANorm += std::norm(v)/(size_t)m_AMC.size();
            INFO("norm(A) = "<<m_ANorm);
            INFO("AMC[0] = "<<m_AMC[0]);

            INFO("Setting up MC for C");
            CoherentSum * C = new CoherentSum(m_SigType.conj(true), m_mps); C->setEvents(mc); C->setMC(mc); C->transferParameters(); C->prepare(); m_CNorm = C->norm(); for (auto& evt:mc) {m_CMC.push_back(C->getValNoCache(evt));} C->reset(true); delete C;

            //for (auto v : m_CMC) m_CNorm += std::norm(v)/(size_t)m_CMC.size();
            INFO("norm(C) = "<<m_CNorm);
            std::vector<complex_t> BMC = {};
            std::vector<complex_t> DMC = {};
            
            INFO("Setting up MC for BD");

            for (int i =0; i < m_TagType.size(); i++){
                INFO("At "<<m_TagType[i]);
                if (!(m_TagType[i]==m_SigType)){
                    
                    INFO("At "<<m_TagType[i]<<" assuming const? "<<m_constAmp[i]);
                    if (m_constAmp[i]){

                        INFO("Assuming constant value for amp so only need 1 event for normalisation");

                        EventList mc_tag =  Generator<>(m_TagType[i], &rndm).generate(1);

                        INFO("Setting up MC for B");
                        CoherentSum * B = new CoherentSum(m_TagType[i].conj(true), m_mps); INFO("size(B) = "<<B->size()); B->setEvents(mc_tag); B->setMC(mc_tag); B->transferParameters(); B->prepare(); m_BNorms.insert(std::pair<size_t, real_t>(i, B->norm())); BMC.push_back(B->getValNoCache(mc_tag[0])); B->reset(true); delete B;

                        
                        for (auto v : BMC) m_BNorms[i] += std::norm(v)/(size_t)BMC.size();
                        INFO("B Norm = "<<m_BNorms[i]);

                        INFO("Setting up MC for D");
                        CoherentSum * D = new CoherentSum(m_TagType[i], m_mps);            D->setEvents(mc_tag); D->setMC(mc_tag); D->transferParameters(); D->prepare(); m_DNorms.insert(std::pair<size_t, real_t>(i, D->norm())); DMC.push_back(D->getValNoCache(mc_tag[0])); D->reset(true); delete D;
                        for (auto v : DMC) m_DNorms[i] += std::norm(v)/(size_t)DMC.size();
                        INFO("D Norm = "<<m_DNorms[i]);
                    }
                    else{
                        INFO("Assuming non-constant value for Amp so generate "<<NInt<<" events for integral");
                        EventList mc_tag =  Generator<>(m_TagType[i], &rndm).generate(NInt);

                        INFO("Setting up MC for B");
                        CoherentSum * B = new CoherentSum(m_TagType[i].conj(true), m_mps); B->setEvents(mc_tag); B->setMC(mc_tag);  B->transferParameters(); B->prepare();  m_BNorms.insert(std::pair<size_t, real_t>(i,B->norm())); for (auto& evt:mc_tag) {BMC.push_back(B->getValNoCache(evt));} B->reset(true); delete B;
                        INFO("Setting up MC for D");
                        CoherentSum * D = new CoherentSum(m_TagType[i], m_mps);            D->setEvents(mc_tag); D->setMC(mc_tag);  D->transferParameters(); D->prepare(); m_DNorms.insert(std::pair<size_t, real_t>(i,D->norm())); for (auto& evt:mc_tag) {DMC.push_back(D->getValNoCache(evt));} D->reset(true); delete D;
                    }
                
                    std::pair<int, std::vector<complex_t> > pB(i, BMC);
                    std::pair<size_t, std::vector<complex_t> > pD(i, DMC);
                    m_BMC.insert(pB);
                    m_DMC.insert(pD);
                }
                //m_DMC.insert(std::pair<EventType, std::vector<complex_t> >(m_TagType[i],DMC));
            }
            //Data
            std::vector<complex_t> vA;
            std::vector<complex_t> vB;
            std::vector<complex_t> vC;
            std::vector<complex_t> vD;
            for (int i=0;i<SigData.size();i++){

                INFO("Setting up data "<<i<<" for A");
                CoherentSum * A = new CoherentSum(m_SigType, m_mps);               A->setEvents(SigData[i]); A->setMC(mc);  A->transferParameters(); A->prepare();for (auto& evt:SigData[i]) {vA.push_back(A->getValNoCache(evt));} A->reset(true); delete A;
                INFO("Setting up data "<<i<<" for B");
                CoherentSum * B = new CoherentSum(m_TagType[i].conj(true), m_mps); B->setEvents(TagData[i]); B->setMC(mc);  B->transferParameters(); B->prepare(); for (auto& evt:TagData[i]) {vB.push_back(B->getValNoCache(evt));} B->reset(true); delete B;
                INFO("Setting up data "<<i<<" for C");
                CoherentSum * C = new CoherentSum(m_SigType.conj(true), m_mps);    C->setEvents(SigData[i]); C->setMC(mc);  C->transferParameters(); C->prepare();for (auto& evt:SigData[i]) {vC.push_back(C->getValNoCache(evt));} C->reset(true); delete C;
                INFO("Setting up data "<<i<<" for D");
                CoherentSum * D = new CoherentSum(m_TagType[i], m_mps);            D->setEvents(TagData[i]); D->setMC(mc);  D->transferParameters(); D->prepare();for (auto& evt:TagData[i]) {vD.push_back(D->getValNoCache(evt));} D->reset(true); delete D;
                PhaseCorrection pc(m_mps);
                INFO("Preparing phase correction for "<<i);
                pc.setEvents(SigData[i]);
                pc.prepareCache();
                m_pc.push_back(pc);
                m_A.push_back(vA);
                m_B.push_back(vB);
                m_C.push_back(vC);
                m_D.push_back(vD);
                if (m_SigType==m_TagType[i]){
                    PhaseCorrection pcT(m_mps);
                    pcT.setEvents(TagData[i]);
                    pcT.prepareCache();
                    m_pcT = pcT;
                    
                }

            }
            for (int i=0; i< GamData.size();i++){
                std::vector<complex_t> vgA;
                std::vector<complex_t> vgC;
                CoherentSum * A = new CoherentSum(m_SigType, m_mps);            A->setEvents(GamData[i]); A->setMC(mc);  A->transferParameters(); A->prepare(); for (auto& evt:GamData[i]) {vgA.push_back(A->getValNoCache(evt));} A->reset(true); delete A;
                CoherentSum * C = new CoherentSum(m_SigType.conj(true), m_mps); C->setEvents(GamData[i]); C->setMC(mc);  C->transferParameters(); C->prepare(); for (auto& evt:GamData[i]) {vgC.push_back(C->getValNoCache(evt));} C->reset(true); delete C;
                PhaseCorrection pcg(m_mps);
                pcg.setEvents(GamData[i]);
                pcg.prepareCache();
                m_pcg.push_back(pcg);
                m_gA.push_back(vgA);
                m_gC.push_back(vgC);
            }


        }
        
        std::vector<complex_t> get_AMC(){ return m_AMC;}
        std::vector<complex_t> get_CMC(){ return m_CMC;}


        real_t NormCorr(int i){
            real_t n=0;
            complex_t w1=0;
            complex_t w2=0;
          //  INFO("At norm");
            for (int j=0; j < m_AMC.size(); j++){
                real_t f = m_pcMC.getValCache(j);
                w1+=m_AMC[j] * std::conj(m_CMC[j]) * exp(complex_t(0,f))/(real_t)m_AMC.size();
                if(m_TagType[i]==m_SigType){
                     w2 += m_CMC[j] * std::conj(m_AMC[j]) * exp(complex_t(0,-f))/(real_t)m_AMC.size();
                }
                else
                {   if(m_constAmp[i]) {
                       w2+=m_BMC[i][0] * std::conj(m_DMC[i][0])/(real_t)m_AMC.size();
                    }
                    else{
                        w2+=m_BMC[i][j] * std::conj(m_DMC[i][j])/(real_t)m_AMC.size();
                    }
                }
            }
           // INFO("w1 = "<<w1);
           // INFO("w2 = "<<w2);
            if (m_TagType[i]==m_SigType) {
                return m_ANorm * m_CNorm * 2 - 2 * std::real(w1* w2);
            }
            else{
                return m_ANorm * m_BNorms[i] + m_CNorm* m_DNorms[i] - 2 * std::real(w1 * w2);
            }
           
             
        }

        real_t LLCorr(int i){
           //     INFO("At LLCorr");
                //auto _LL =  make_likelihood( m_SigData[i], m_TagData[i], m_Psi[i]
//                pc.stop();
            real_t ll =0;
            
            real_t n = NormCorr(i);
        //    INFO("norm = "<<n);
            for (size_t j=0; j < m_A.size(); j++){
                real_t f1 = m_pc[i].getValCache(j);
                if (m_SigType==m_TagType[i]){
                    real_t f2 = m_pcT.getValCache(j);
                    real_t pdf = std::norm(m_A[i][j] * m_B[i][j]) + std::norm(m_C[i][j] * m_D[i][j]) - 2 * std::real(m_A[i][j] * m_B[i][j] * std::conj(m_C[i][j] * m_D[i][j])  *exp(complex_t(0, f1 - f2)) );
                    ll+=log(pdf/n);
                }
                else{
                    real_t pdf = std::norm(m_A[i][j] * m_B[i][j]) + std::norm(m_C[i][j] * m_D[i][j]) - 2 * std::real(m_A[i][j] * m_B[i][j] * std::conj(m_C[i][j] * m_D[i][j])  *exp(complex_t(0, f1)) );
                    ll+=log(pdf/n);

                }
            }
            

//            if (m_debug) INFO("LL "<<i<<" = "<<_LL<<" took "<<pc.t_duration<<" for "<<m_SigData[i].size()<<" events");
            return -2*ll;

        }

         complex_t sumFactor(int gamSign, bool useXY){
            if (useXY==1){
                if (gamSign > 0) return complex_t(m_mps["CKM::x+"]->mean(), m_mps["CKM::y+"]->mean());
                if (gamSign < 0) return complex_t(m_mps["CKM::x-"]->mean(), m_mps["CKM::y-"]->mean());
            }
            else{
                return (complex_t)m_mps["CKM::rB"]->mean() * exp(complex_t(0, m_mps["CKM::dB"]->mean() + gamSign * m_mps["CKM::gamma"]->mean()));
            }
          }

        real_t NormGam(size_t i){
            real_t n=0;
            complex_t sf = sumFactor(m_gammaSigns[i], m_useXYs[i]);
            complex_t w1=0;
            complex_t w2=0;
            for (size_t j=0; j<m_AMC.size(); j++){
                real_t f = m_pcMC.getValCache(j);
                w1 += m_AMC[j] * std::conj(m_CMC[j] * sf) * cos(f)/(real_t)m_AMC.size(); 
                w2 += -m_AMC[j] * std::conj(m_CMC[j] * sf) * sin(f)/(real_t)m_AMC.size(); 
            }
            if (m_BConj[i]==1){
                return m_CNorm + std::norm(sf) * m_ANorm + 2 * std::real( std::conj(w1) + complex_t(0,1) *  std::conj(w2) );
            }
            else{
                return m_ANorm + std::norm(sf) * m_CNorm + 2 * std::real(w1 + complex_t(0,1) * w2);
            }
        }

        real_t LLGam(size_t i){
            real_t ll=0;
//            INFO("At LL Gam");
            real_t n = NormGam(i);
 //           INFO("norm = "<<n);
            complex_t sf = sumFactor(m_gammaSigns[i], m_useXYs[i]);
            for (size_t j=0;j<m_gA[i].size(); j++){
                real_t pdf = 0;

                real_t f = m_pcg[i].getValCache(j);
               
                if (m_BConj[i]==1){

                    real_t pdf = std::norm(m_gC[i][j]) + std::norm(sf * m_gA[i][j]) + 2 * std::real( m_gC[i][j] * std::conj(sf * m_gA[i][j]) * exp(complex_t(0, f))  );
                    ll+=log(pdf/n);
                }
                else{
                    real_t pdf = std::norm(m_gA[i][j]) + std::norm(sf * m_gC[i][j]) + 2 * std::real( m_gA[i][j] * std::conj(sf * m_gC[i][j]) * exp(complex_t(0, f))  );
                    ll+=log(pdf/n);
                }
            }
            return -2 * ll;
        }

        real_t LL(){
            real_t ll=0;
            for (size_t i=0;i<m_gA.size();i++){
                ll+=LLGam(i);
            }
            for (size_t i=0;i<m_A.size();i++){
                ll+=LLCorr(i);
            }
            return ll;
        }

        MinuitParameterSet getMPS(){
            return m_mps;
        }



        double getVal(){
            double ll =0 ;
            for (int i=0;i<m_A.size();i++){
                ll+=LLCorr(i);
            }
            for (int i=0; i< m_gA.size();i++){
                ll+=LLGam(i);
            }
            return ll;

           
        }

  };
}
#endif
