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
        std::vector<EventList> m_SigData;
        std::vector<EventList> m_TagData;
        std::vector<EventList> m_GamData;


        std::vector<EventList> m_SigInt;
        std::vector<EventList> m_TagInt;
        std::vector<EventList> m_GamInt;

        std::vector<EventType> m_SigType;
        std::vector<EventType> m_TagType;
        std::vector<EventType> m_GamType;

        MinuitParameterSet m_mps;
        bool m_debug;
        std::vector<pCorrelatedSum> m_Psi;
        std::vector<pCoherentSum> m_A;
        std::vector<std::string> m_SumFactors;
        std::vector<int> m_gammaSigns;
        std::vector<bool> m_useXYs;
        std::vector<bool> m_BConj;
        


      public:
        CombGamCorrLL() = default;
        CombGamCorrLL(std::vector<EventList> SigData, 
               std::vector<EventList> TagData, 
               std::vector<EventList> GamData, 

                   std::vector<EventList> SigInt, 
                   std::vector<EventList> TagInt, 
                   std::vector<EventList> GamInt, 
                   

                   std::vector<EventType> SigType,
                   std::vector<EventType> TagType,
                   std::vector<EventType> GamType,

                   MinuitParameterSet mps,
                   std::vector<std::string> sumFactors,
                   std::vector<int> gammaSigns,
                   std::vector<bool> useXYs,
                   std::vector<bool> BConj):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        m_GamData(GamData),

                        m_SigInt(SigInt),
                        m_TagInt(TagInt),
                        m_GamInt(GamInt),

                        m_SigType(SigType),
                        m_TagType(TagType),
                        m_GamType(GamType),

                        m_mps(mps),
                        m_SumFactors(sumFactors),
                        m_gammaSigns(gammaSigns),
                        m_useXYs(useXYs),
                        m_BConj(BConj),
                        m_debug(NamedParameter<bool>("CombGamCorrLL::Debug", false, "Debug CombLL"))
                        {
                            std::vector<pCoherentSum> pCS = {};
                            for (auto i=0; i < m_GamData.size() ; i++){
                                pCoherentSum _pCS = pCoherentSum(m_GamType[i], m_mps, m_SumFactors[i], m_gammaSigns[i], m_useXYs[i], m_BConj[i]);
                                _pCS.setEvents(m_GamData[i]);
                                _pCS.setMC(m_GamInt[i]);
                                _pCS.prepare();
                                m_A.push_back(_pCS);
                            }
                            for (auto i=0; i < m_SigData.size() ; i++){
                                pCorrelatedSum _pCS =  pCorrelatedSum(m_SigType[i], m_TagType[i], m_mps);


                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                _pCS.setMC(m_SigInt[i], m_TagInt[i]);
                                _pCS.prepare();
                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        double LL_Corr(int i){

                //auto _LL =  make_likelihood( m_SigData[i], m_TagData[i], m_Psi[i]);
                auto psi = m_Psi[i];
                auto sigDat = m_SigData[i];
                auto tagDat = m_TagData[i];
                ProfileClock pcN;
                pcN.start();
                auto norm = psi.norm();
                pcN.stop();
                if (m_debug) INFO("LL "<<i<<" norm = "<<norm<<" took "<<pcN.t_duration<<" for "<<m_SigInt[i].size()<<" integration events");
                real_t _LL = 0;
                ProfileClock pc;
                pc.start();

                #pragma omp parallel for reduction( +: _LL )
                for (int j=0;j<sigDat.size();j++){
                    auto evt1 = sigDat[j];
                    auto evt2 = tagDat[j];
                    auto prob = std::norm(psi.getVal(evt1,evt2))/norm;
                    _LL += -2 * log(prob);


                }
                pc.stop();

            if (m_debug) INFO("LL "<<i<<" = "<<_LL<<" took "<<pc.t_duration<<" for "<<m_SigData[i].size()<<" events");
            return _LL;

        }
        MinuitParameterSet getMPS(){
            return m_mps;
        }

        double LL_Gam(int i){
            auto psi = m_A[i];
            auto sigDat = m_GamData[i];
            ProfileClock pcN;
            pcN.start();
            auto norm = psi.norm();
            pcN.stop();
            if (m_debug) INFO("LL "<<i<<" norm = "<<norm<<" took "<<pcN.t_duration<<" for "<<m_GamInt[i].size()<<" integration events");
            real_t _LL=0;
            ProfileClock pc;
            pc.start();
            #pragma omp parallel for reduction( +: _LL )
            for (int j=0;j<sigDat.size();j++){
                    auto evt1 = sigDat[j];
                    auto prob = std::norm(psi.getVal(evt1))/norm;
                    _LL += -2 * log(prob);
            }
            pc.stop();
            if (m_debug) INFO("LL "<<i<<" = "<<_LL<<" took "<<pc.t_duration<<" for "<<m_GamData[i].size()<<" events");
            

            return _LL;

        }


        double getVal(){
            double ll =0 ;

            for (auto i=0; i < m_SigData.size(); i++){
                if (m_debug) INFO("Calculating "<<i+1<<"/"<<m_SigData.size()<<" LL corr");
                ll += LL_Corr(i);
            }
            
            for (auto i=0; i < m_GamData.size(); i++){
              
                if (m_debug) INFO("Calculating "<<i+1<<"/"<<m_GamData.size()<<" LL gamma");

                ll += LL_Gam(i);
            }
            if (m_debug) INFO("LL = "<<ll);
            return ll;
        }

  };
}
#endif
