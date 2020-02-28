#ifndef SUM_LL
#define SUM_LL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/CorrelatedLL.h"

#include <tuple>

namespace AmpGen
{
  template <class LL>
  class SumLL
    {
    private : 
        std::vector<LL> m_LLs;
        bool m_debug;
    public:
        //Default Constructor
        SumLL() = default;

        //Build from a std::vector of LogLikelihoods
        SumLL(std::vector<LL> LLs) : 
            m_LLs(LLs),
            m_debug(NamedParameter<bool>("SumLL::debug", false, "Flag to debug SumLL"))
        {}

        double getVal(){
            double res=0;
            for(auto ll : m_LLs){
                if (m_debug){
                    INFO("LL = "<<ll.getVal());
                }
                res += ll.getVal();
            }
            return res;
        }
    };


       
  
  
}
#endif
