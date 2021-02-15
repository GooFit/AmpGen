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
      polyLASSO(SimFit ll, const MinuitParameterSet mps) : m_mps(mps), m_ll(ll), m_order(NamedParameter<size_t>( "pCorrelatedSum::Order") ), m_debug(NamedParameter<bool>("LASSO::debug", false)), m_start(NamedParameter<size_t>("pCorrelatedSum::Start", 0)) {}
      double getVal()
      {
        double LL = m_ll.getVal();
        double pen = penalty();
       
        return LL + pen;
      }

      real_t penaltyPerOrder(size_t order){
        real_t penalty=0;
        for (size_t i=0;i<order+1;i++){
            size_t i1 = i;
            size_t i2 = order - i;
            penalty += std::abs(m_mps["pCorrelatedSum::C" + std::to_string(i1) + std::to_string(i2)]->mean());
        }
        return penalty;
      }

      double penalty(){
        real_t lambda = m_mps["LASSO::lambda"]->mean() ;
        if (m_debug) INFO("lambda = "<<lambda);
        double pen = 0;
        for (size_t i=m_start;i<m_order;i++){
          pen += penaltyPerOrder(i);
        }
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