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
      polyLASSO(SimFit ll, const MinuitParameterSet mps) : m_mps(mps), m_ll(ll), m_order(NamedParameter<size_t>( "pCorrelatedSum::Order") ), m_lambda(NamedParameter<real_t>("LASSO::Lambda")), m_debug(NamedParameter<bool>("LASSO::debug", false)) {}
      double getVal()
      {
        double LL = m_ll.getVal();
        double pen = penalty();
       
        return LL + pen;
      }

      double penalty(){

        if (m_debug) INFO("lambda = "<<m_lambda);

        double pen = 0;
        for (size_t i=0;i<m_order;i++){
            for (size_t j=0; j<m_order - i; j++){
                pen += std::abs(m_mps["pCorrelatedSum::C" + std::to_string(i) + std::to_string(j)]->mean());
            }
        }
        return m_lambda * pen;

      }

   
    private:
        SimFit m_ll;
        MinuitParameterSet m_mps;
        size_t m_order;
        real_t m_lambda;
        bool m_debug;


  };
} // namespace AmpGen

#endif