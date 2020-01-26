#ifndef AMPGEN_POLARISEDAMPLITUDE_H
#define AMPGEN_POLARISEDAMPLITUDE_H

#include <stddef.h>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <array>

#include "AmpGen/Types.h"
#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator2.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/MinuitParameter.h"

#include "TMatrixD.h"

namespace AmpGen
{
  class LinearErrorPropagator;
  class MinuitParameterSet;
  class FitFraction;
  class Event;
  class EventList;
  class MinuitProxy; 

  class PolarisedSum
  {
    public: 
      typedef Integrator<10>        integrator;

      PolarisedSum() = default; 
      PolarisedSum(const EventType&, AmpGen::MinuitParameterSet&, const std::vector<MinuitProxy>& = {});
      void prepare();
      void setEvents(AmpGen::EventList&);
      void setMC(AmpGen::EventList&);
      void reset(const bool& = false);
      void debug(const AmpGen::Event&);
      void debug_norm(); 
      void setWeight(MinuitProxy);
      double getWeight() const;
      void calculateNorms(const std::vector<bool>&); 
      void generateSourceCode(const std::string&, const double& = 1, bool = false);
      void build_probunnormalised();
      Expression probExpression(const Tensor&, const std::vector<Expression>&, DebugSymbols* = nullptr) const; 
      size_t size() const;  
      real_t norm() const;
      complex_t norm(const size_t&, const size_t&, integrator* = nullptr); 
      inline real_t operator()(const AmpGen::Event& evt) const { return m_weight * prob_unnormalised(evt) / m_norm; }
      real_t prob_unnormalised(const AmpGen::Event&) const;
      real_t prob(const AmpGen::Event&) const;
      real_t getValNoCache(const AmpGen::Event&) ;
      std::vector<FitFraction> fitFractions(const LinearErrorPropagator&);
      std::vector<TransitionMatrix<std::vector<complex_t>>> matrixElements() const;
      void transferParameters(); 
      Tensor transitionMatrix();
      TransitionMatrix<std::vector<complex_t>> operator[](const size_t& i) const { return m_matrixElements[i] ; } 

    private: 
      size_t                        m_nCalls      = {0};
      real_t                        m_norm        = {1};
      EventList*                    m_events      = {nullptr};
      MinuitParameterSet*           m_mps         = {nullptr};
      MinuitProxy                   m_weight      = {nullptr,1}; 
      std::vector<MinuitProxy>      m_pVector     = {}; 
      bool                          m_verbosity   = {0};
      bool                          m_debug       = {0};
      integrator                    m_integrator;
      std::vector<Bilinears>        m_norms;
      EventType                     m_eventType;
      std::string                   m_prefix      = "";
      std::vector<complex_t>        m_rho;
      std::vector<size_t>           m_integIndex; 
      AmplitudeRules                m_rules;  
      std::pair<unsigned, unsigned> m_dim; 
      std::vector<TransitionMatrix<std::vector<complex_t>>>       m_matrixElements;  
      CompiledExpression<real_t, const real_t*, const complex_t*> m_probExpression; 

      std::vector<std::vector<int>> indexProduct(const std::vector<std::vector<int>>&, const std::vector<int>&) const;
      std::vector<int> polarisations(const std::string&) const ;
  };
} // namespace AmpGen

#endif
