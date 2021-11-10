#ifndef POLYLASSO_H_
#define POLYLASSO_H_
#include "AmpGen/SimPDF.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
namespace AmpGen
{
  class polyLASSO
  {
    public:
      polyLASSO(SimFit ll, const MinuitParameterSet mps) : m_mps(mps), m_ll(ll), m_order(NamedParameter<size_t>( "PhaseCorrection::Order") ), m_debug(NamedParameter<bool>("LASSO::debug", false)), m_start(NamedParameter<size_t>("PhaseCorrection::Start", 0)) {}
      double getVal()
      {
        double LL = m_ll.getVal();
        double pen = penalty();
       
        return LL + pen;
      }

      real_t penaltyPhaseCorrection(){
        real_t penalty=0;
        for (size_t i=0;i<m_order+1;i++){
            for (size_t j=0;j<m_order+1-i;j++){
                int i1 = i;
                int i2 = 2 * j + 1;
                penalty += std::abs(m_mps["PhaseCorrection::C" + std::to_string(i1) + "_" + std::to_string(i2)]->mean());
            }
        }
        return penalty;
      }

      double penalty(){
        real_t lambda = m_mps["LASSO::lambda"]->mean() ;
        if (m_debug) INFO("lambda = "<<lambda);
        double pen = penaltyPhaseCorrection();
//        for (auto& p : m_mps){
//          if (p->isFree()) pen += std::abs(p->mean());
//        }
       
        return  lambda * pen;

      }

   
    private:
        SimFit m_ll;
        MinuitParameterSet m_mps;
        size_t m_order;
        size_t m_start;
        bool m_debug;


  };
} // namespace AmpGen

#endif