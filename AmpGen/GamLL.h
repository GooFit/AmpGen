#ifndef GAMLL_H_
#define GAMLL_H_

#include "AmpGen/CoherentSum.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Generator.h"
#include "TRandom3.h"
#include "AmpGen/MinuitParameterSet.h"

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"

namespace AmpGen
{
    class GamLL
  {
      private:
        std::vector<std::vector<complex_t > > m_A = {};
        std::vector<std::vector<complex_t > > m_C = {};
        std::vector<complex_t> m_AMC = {};
        std::vector<complex_t> m_CMC = {};
        std::vector<int> m_gamSign;
        std::vector<int> m_useXY;
        std::vector<int> m_conj;
        real_t m_normA = 0;
        real_t m_normC = 0;

        MinuitParameterSet m_mps;
        PhaseCorrection m_pcMC;
        std::vector<PhaseCorrection> m_pc = {};
      
      public:
        complex_t sumFactor(int gamSign, bool useXY){
            if (useXY==1){
                if (gamSign > 0) return complex_t(m_mps["CKM::x+"]->mean(), m_mps["CKM::y+"]->mean());
                if (gamSign < 0) return complex_t(m_mps["CKM::x-"]->mean(), m_mps["CKM::y-"]->mean());
            }
            else{
                return (complex_t)m_mps["CKM::rB"]->mean() * exp(complex_t(0, m_mps["CKM::dB"]->mean() + gamSign * m_mps["CKM::gamma"]->mean()));
            }
          }
        
        std::vector<complex_t> get_AMC() { return m_AMC;}
        std::vector<complex_t> get_CMC() { return m_CMC;}

        std::vector<std::vector<complex_t> > get_A() { return m_A; }
        std::vector<std::vector<complex_t> > get_C() { return m_C; }
        real_t get_normA() { return m_normA; }
        real_t get_normC() { return m_normC; }
        real_t norm(int i) {
            real_t nA = m_normA;
            real_t nC = m_normC;
            std::vector<complex_t> w = {0,0};


            for (int j=0; j<m_AMC.size(); j++){
                real_t f = m_pcMC.getValCache(j);
                complex_t r = m_AMC[j] * std::conj(m_CMC[j])/(real_t)m_AMC.size();
                w[0]+=r * cos(f);
                w[1]+=-r * sin(f);
            }
            complex_t sf = sumFactor(m_gamSign[i], m_useXY[i]);
            real_t inter = 2 * std::real( std::conj(sf) * (w[0] + complex_t(0,1) * w[1]) );
            if (m_conj[i]==1){
                inter = 2 * std::real(std::conj(sf) * (std::conj(w[0]) + complex_t(0,1) * std::conj(w[1]) ) );
                return nC + std::norm(sf) * nA+ inter;
            }
            return nA + std::norm(sf) * nC+ inter;
        }

        real_t LL(int i){
            real_t ll =0;
            real_t n = norm(i);
            for (int j=0; j < m_A[i].size(); j++){

                real_t f = m_pc[i].getValCache(j);

                if (m_conj[i]==1){

                    ll += log(std::abs(m_C[i][j]* exp(complex_t(0, f/2)) + sumFactor(m_gamSign[i], m_useXY[i]) * m_A[i][j] * exp(complex_t(0, -f/2)))/n);
                }
                else{
                    
                    ll += log(std::abs(m_A[i][j] * exp(complex_t(0, f/2)) + sumFactor(m_gamSign[i], m_useXY[i]) * m_C[i][j] * exp(complex_t(0, -f/2)))/n);
                    
                }

            }
            return -2 * ll;
        }

        real_t operator()(){
            real_t ll=0;
            for (int i =0; i<m_A.size(); i++){
                ll+=LL(i);
            }
            return ll;
        }
        

        GamLL(std::vector<EventList> SigData, 
              EventType SigType,
              MinuitParameterSet mps,
              std::vector<int> gammaSigns,
              std::vector<int> useXYs,
              std::vector<int> conjs
             ):
  
              m_gamSign(gammaSigns),
              m_useXY(useXYs),
              m_conj(conjs),
              m_mps(mps),
              m_pcMC(PhaseCorrection(mps))
                
                {
                    TRandom3 rndm;
                    rndm.SetSeed( NamedParameter<size_t>("GamLL::Seed", 0) );
                    gRandom = &rndm;
                    size_t NInt = NamedParameter<size_t>("GamLL::NInt", 1e5);
                    EventList mc =  Generator<>(SigType, &rndm).generate(NInt);
                    m_pcMC.setEvents(mc);
                    m_pcMC.prepareCache();
                    for (int i=0; i < SigData.size(); i++){

                        PhaseCorrection pc(mps);
                        pc.setEvents(SigData[i]);
                        pc.prepareCache();
                        m_pc.push_back(pc);


                        std::vector<complex_t> AVals_i;
                        CoherentSum *  A = new CoherentSum(SigType, mps);
                        A->setEvents(SigData[i]);
                        A->setMC(mc);
                        A->transferParameters();
                        A->prepare();

                        for (auto& evt : SigData[i]){
                            AVals_i.push_back(A->getVal(evt));
                        }
                        m_A.push_back(AVals_i);
                        
                        A->reset(true);
                        delete A;

                        std::vector<complex_t> CVals_i;
                        CoherentSum *  C = new CoherentSum(SigType.conj(true), mps);
                        C->setEvents(SigData[i]);
                        C->setMC(mc);
                        C->transferParameters();
                        C->prepare();

                        for (auto& evt : SigData[i]){
                            CVals_i.push_back(C->getVal(evt));
                        }
                        m_C.push_back(CVals_i);
                        
                        C->reset(true);
                        delete C;
                    }

                    CoherentSum *  A = new CoherentSum(SigType, mps);
                    A->setEvents(mc);
                    A->setMC(mc);
                    A->transferParameters();
                    A->prepare();
                    m_normA = A->norm();

                    for (auto& evt : mc){
                        m_AMC.push_back(A->getVal(evt));
                    }
                    
                    A->reset(true);
                    delete A;


                    CoherentSum *  C = new CoherentSum(SigType.conj(true), mps);
                    C->setEvents(mc);
                    C->setMC(mc);
                    C->transferParameters();
                    C->prepare();

                    m_normC = C->norm();

                    for (auto& evt : mc){
                        m_CMC.push_back(C->getVal(evt));
                    }
                    
                    C->reset(true);
                    delete C;
                }

  };
       

}

#endif