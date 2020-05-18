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
        


      public:
        CombCorrLL() = default;
        CombCorrLL(std::vector<EventList> SigData, 
                   std::vector<EventList> TagData, 
                   std::vector<EventList> SigInt, 
                   std::vector<EventList> TagInt,
                   std::vector<EventType> SigType,
                   std::vector<EventType> TagType,
                   MinuitParameterSet mps):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        m_SigInt(SigInt),
                        m_TagInt(TagInt),
                        m_SigType(SigType),
                        m_TagType(TagType),
                        m_mps(mps),
                        m_debug(NamedParameter<bool>("CombCorrLL::Debug", false, "Debug CombCorrLL"))
                        {
                            std::vector<pCorrelatedSum> pCS = {};
                            for (auto i=0; i < m_SigData.size() ; i++){
                                pCorrelatedSum _pCS = pCorrelatedSum(m_SigType[i], m_TagType[i], m_mps);
                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                _pCS.setMC(m_SigInt[i], m_TagInt[i]);
                                _pCS.prepare();
                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        double LL(int i){
            auto _LL =  make_likelihood( m_SigData[i], m_TagData[i] , m_Psi[i]);
            return _LL.getVal();
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
