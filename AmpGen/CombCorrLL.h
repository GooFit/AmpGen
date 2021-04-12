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
        std::vector<std::string> m_SumFactors;
        


      public:
        CombCorrLL() = default;
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
                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        double LL(int i){

                //auto _LL =  make_likelihood( m_SigData[i], m_TagData[i] , m_Psi[i]);
                auto psi = m_Psi[i];
                real_t norm = psi.norm();
                real_t ll = 0;
                #pragma omp parallel for reduction( +: ll )
                for ( unsigned int j = 0; j < m_SigData[i].size(); ++j ) {
                  auto evt1 = m_SigData[i][j];
                  auto evt2 = m_TagData[i][j];
                  auto prob = std::norm(psi.getVal(evt1,evt2))/norm;
                  ll += -2* log(prob);

                }
            if (m_debug) INFO("LL = "<<ll);
            if (m_debug){
              auto norm = psi.norm();
              INFO("norm = "<<norm);
            }
            return ll;

        }

        double getVal(){
            double ll =0 ;

            #pragma omp parallel for reduction( +: ll )
            for (auto i=0; i < m_SigData.size(); i++){
                if (m_debug) INFO("LL_"<<i<<" = "<<LL(i));
                ll += LL(i);
            }
            if (m_debug) INFO("LL = "<<ll);
            return ll;
        }

  };
}
#endif
