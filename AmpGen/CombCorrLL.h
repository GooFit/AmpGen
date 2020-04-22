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
                        m_mps(mps)
                        {

                        }
        double LL(int i){
            pCorrelatedSum pCS = pCorrelatedSum(m_SigType[i], m_TagType[i], m_mps);
            pCS.setEvents(m_SigData[i], m_TagData[i]);
            pCS.setMC(m_SigInt[i], m_TagInt[i]);
            pCS.prepare();
            double ll =0 ;
            for (int j=0; j<m_SigData.size();j++){
                double Psi2 = pCS.prob(m_SigData[i][j], m_TagData[i][j]);
                ll += log(Psi2);
            }
            return ll;
        }

        double getVal(){
            double ll =0 ;
            for (int i=0; i < m_SigData.size(); i++){
                ll += LL(i);
            }
            return -2 * ll;
        }

        double operator() (){
            return getVal();
        }

  };
}
#endif