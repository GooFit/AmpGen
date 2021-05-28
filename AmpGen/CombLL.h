#ifndef COMBLL
#define COMBLL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/SumPDF.h"

#include <tuple>
namespace AmpGen
{
  class EventList;
  /**
   * @class CombLL
   * @brief A combined log-likelihood for elated amplitudes.
   **/
  class CombLL
  {
      private:
        std::vector<EventList> m_SigData;

       
        EventList m_SigInt;

        std::vector<EventType> m_SigType;

        MinuitParameterSet m_mps;
        bool m_debug;
        std::vector<pCoherentSum> m_Psi;
        std::vector<std::string> m_SumFactors;
        std::vector<int> m_gammaSigns;
        std::vector<int> m_useXYs;
        std::vector<int> m_conj;
        


      public:
        CombLL() = default;
        CombLL(std::vector<EventList> SigData, 

                   EventList SigInt, 

                   std::vector<EventType> SigType,

                   MinuitParameterSet mps,
                   std::vector<std::string> sumFactors,
                   std::vector<int> gammaSigns,
                   std::vector<int> useXYs,
                   std::vector<int> conjs
                   ):
                        m_SigData(SigData),

                        m_SigInt(SigInt),

                        m_SigType(SigType),

                        m_mps(mps),
                        m_SumFactors(sumFactors),
                        m_gammaSigns(gammaSigns),
                        m_useXYs(useXYs),
                        m_conj(conjs),
                        m_debug(NamedParameter<bool>("CombLL::Debug", false, "Debug CombLL"))
                        {
                            std::vector<pCoherentSum> pCS = {};
                            for (auto i=0; i < m_SigData.size() ; i++){
                                pCoherentSum _pCS = pCoherentSum(m_SigType[i], m_mps, m_gammaSigns[i], m_useXYs[i], m_conj[i]);
                                _pCS.setEvents(m_SigData[i]);
                                _pCS.setMC(m_SigInt);
                                _pCS.prepare();
                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        double LL(int i){
            auto psi = m_Psi[i];
            auto sigDat = m_SigData[i];
            ProfileClock pcN;
            
            pcN.start();
            auto norm = psi.norm();
            pcN.stop();
            if (m_debug) INFO("LL "<<i<<" norm = "<<norm<<" took "<<pcN.t_duration<<" for "<<m_SigInt[i].size()<<" integration events");
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
            if (m_debug) INFO("LL "<<i<<" = "<<_LL<<" took "<<pc.t_duration<<" for "<<m_SigData[i].size()<<" events");
            

            return _LL;






        }

        double getVal(){
            double ll =0 ;

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
