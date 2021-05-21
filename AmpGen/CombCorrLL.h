#ifndef COMBCORRLL
#define COMBCORRLL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/Generator.h"

#include <tuple>
namespace AmpGen
{
  class EventList;
  /**
   * @class CombCorrLL
   * @brief A combined log-likelihood for correlated amplitudes.
   **/
  class CombCorrLL
  {
      private:
        std::vector<EventList> m_SigData;
        std::vector<EventList> m_TagData;
        std::vector<EventList> m_SigInt;
        std::vector<EventList> m_TagInt;
        std::vector<EventType> m_SigType;
        std::vector<EventType> m_TagType;
        MinuitParameterSet m_mps;
        bool m_debug;
        std::vector<pCorrelatedSum> m_Psi;
        std::vector<pCorrelatedSum*> m_Psi2;
        std::vector<std::string> m_SumFactors;
        


      public:
        CombCorrLL() = default;
        CombCorrLL(std::vector<pCorrelatedSum>psi):
          m_Psi(psi)
          {
            for(int i=0;i<psi.size();i++){
              INFO("psi["<<i<<"]  ");//<<psi[i]);

             
              
            }

          }

        void setEvents(int i, EventList list1, EventList list2){
          m_Psi[i].setEvents(list1, list2);
        }
        void setMC(int i, EventList list1, EventList list2){
          m_Psi[i].setMC(list1, list2);
        }

        

        CombCorrLL(std::vector<EventList> SigData, 
                   std::vector<EventList> TagData, 
               
               
                   EventType SigType,
                   std::vector<EventType> TagType,
                   MinuitParameterSet mps,
                   size_t NInt,
                   size_t seed
        ):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        //m_SigInt(SigInt),
               
                        //m_SigType(SigType),
                        m_TagType(TagType),
                        m_mps(mps),
                       // m_SumFactors(sumFactors),
                        m_debug(NamedParameter<bool>("CombCorrLL::Debug", false, "Debug CombCorrLL"))
                        {
          TRandom3 rndm;
          rndm.SetSeed( seed );
          gRandom = &rndm;


                           
                            m_SigType.push_back(SigType);

                            auto mcSig = Generator<>(SigType, &rndm).generate(NInt);

                            std::vector<pCorrelatedSum> pCS = {};
                            for (auto i=0; i < m_SigData.size() ; i++){
                              if (m_debug) INFO("Making pCorrelatedSum");
                                pCorrelatedSum _pCS = pCorrelatedSum(SigType, m_TagType[i], m_mps);

                                auto mcTag = Generator<>(TagType[i], &rndm).generate(NInt);
                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                _pCS.setMC(mcSig, mcTag);
                                if (m_debug) INFO("Preparing "<<i+1<<"/"<<m_SigData.size());
                                _pCS.prepare();
                                if (m_debug) INFO("Prepared "<<i+1<<"/"<<m_SigData.size());
                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
     

        CombCorrLL(std::vector<EventList> SigData, 
                   std::vector<EventList> TagData, 
                   std::vector<EventList> SigInt, 
                   std::vector<EventList> TagInt,
                   std::vector<EventType> SigType,
                   std::vector<EventType> TagType,
                   MinuitParameterSet mps,
                   std::vector<std::string> sumFactors):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        m_SigInt(SigInt),
                        m_TagInt(TagInt),
                        m_SigType(SigType),
                        m_TagType(TagType),
                        m_mps(mps),
                        m_SumFactors(sumFactors),
                        m_debug(NamedParameter<bool>("CombCorrLL::Debug", false, "Debug CombCorrLL"))
                        {
                            std::vector<pCorrelatedSum> pCS = {};
                            for (auto i=0; i < m_SigData.size() ; i++){
                                pCorrelatedSum _pCS = pCorrelatedSum(m_SigType[i], m_TagType[i], m_mps, m_SumFactors[i]);
                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                _pCS.setMC(m_SigInt[i], m_TagInt[i]);
                                if (m_debug) INFO("Preparing "<<i+1<<"/"<<m_SigData.size());
                                _pCS.prepare();
                                if (m_debug) INFO("Prepared "<<i+1<<"/"<<m_SigData.size());
                                if (m_debug) INFO("Norm = "<<_pCS.norm());
                                if (m_debug) INFO("LL = "<<_pCS.LL());

                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        double LL(int i){

                //auto _LL =  make_likelihood( m_SigData[i], m_TagData[i] , m_Psi[i]);
                auto psi = m_Psi[i];
                if (m_debug) INFO("Getting norm for "<<i);
                real_t norm = psi.norm();
                if (m_debug) INFO("LL = "<<psi.LL());
                if (m_debug) INFO("Norm = "<<norm);

                real_t ll = 0;
                #pragma omp parallel for reduction( +: ll )
                for ( unsigned int j = 0; j < m_SigData[i].size(); ++j ) {
                  auto evt1 = m_SigData[i][j];
                  auto evt2 = m_TagData[i][j];
                  auto prob = std::norm(psi.getVal(evt1,evt2))/norm;
                  ll += -2* log(prob);

                }
            if (m_debug) INFO("LL = "<<ll);

            return ll;

        }

        double getVal(){
          real_t ll=0;
          for (int i=0;i<m_Psi.size();i++) ll+=LL(i);

          return ll;

        }

  };
}
#endif
