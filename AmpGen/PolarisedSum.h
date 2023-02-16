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
#include "AmpGen/Integrator.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/MinuitParameter.h"

namespace AmpGen
{
  class LinearErrorPropagator;
  class MinuitParameterSet;
  class FitFraction;
  class MinuitProxy; 

  class PolarisedSum
  {
    public: 
    #if ENABLE_AVX
      using EventList_type  = EventListSIMD; 
    #else 
      using EventList_type  = EventList; 
    #endif

      PolarisedSum() = default; 
      PolarisedSum(const EventType&, MinuitParameterSet&, const std::vector<MinuitProxy>& = {});
      void prepare();
      void setEvents(EventList_type&);
      void setMC(EventList_type&);
      #if ENABLE_AVX
        void setEvents(EventList& evts){ m_ownEvents = true; setEvents( *new EventList_type(evts)) ; };
        void setMC(EventList& evts){ setMC( *new EventList_type(evts)) ; };
        double operator()(const double*, const unsigned) const; 
      #endif
      real_v operator()(const real_v*, const unsigned) const; 
      real_t  operator()(const Event& evt) const;
      void reset(const bool& = false);
      void debug(const Event&);
      void debug_norm(); 
      void setWeight(MinuitProxy);
      double getWeight() const;
      void updateNorms(); 
      void generateSourceCode(const std::string&, const double& = 1, bool = false);
      Expression probExpression(const Tensor&, const std::vector<Expression>&, const std::vector<Expression>&, DebugSymbols* = nullptr) const; 
      size_t size() const;  
      real_t norm() const;
      complex_t norm(const size_t&, const size_t&, Integrator* = nullptr); 
      real_t getValNoCache(const Event&) const;
      std::vector<FitFraction> fitFractions(const LinearErrorPropagator&);
      std::vector<MatrixElement> matrixElements() const;
      void transferParameters(); 
      Tensor transitionMatrix() const;
      const MatrixElement& operator[](const size_t& i) const { return m_matrixElements[i] ; } 
      std::function<real_t(const Event&)> evaluator(const EventList_type* = nullptr) const; 
      KeyedFunctors<double(Event)> componentEvaluator(const EventList_type* = nullptr) const;     
      EventType eventType() const{ return m_eventType; }
      ~PolarisedSum();
    private: 
      size_t                        m_nCalls      = {0};
      real_t                        m_norm        = {1};
      EventList_type*               m_events      = {nullptr};
      FunctionCache<EventList_type, complex_v, Alignment::AoS> m_cache    = {};
      FunctionCache<FunctionCache<EventList_type, complex_v, Alignment::AoS>, real_v, Alignment::AoS>  m_pdfCache = {}; 
      bool                          m_ownEvents   = {false};
      MinuitParameterSet*           m_mps         = {nullptr};
      MinuitProxy                   m_weight      = {nullptr,1}; 
      MinuitParameter*              m_polParam    = {nullptr};
      std::vector<MinuitProxy>      m_pVector     = {};
      std::vector<MinuitProxy>      m_pfVector = {};  
      bool                          m_verbosity   = {0};
      bool                          m_debug       = {0};
      Integrator                    m_integrator;
      std::vector<Bilinears>        m_norms;
      EventType                     m_eventType;
      std::string                   m_prefix      = "";
      std::vector<complex_t>        m_rho;
      std::vector<size_t>           m_integIndex; 
      std::pair<unsigned, unsigned> m_dim; 
      std::vector<MatrixElement>                                   m_matrixElements;  
//      CompiledExpression<void(real_t*, const std::size_t*, const real_t*, const complex_v*)>  m_probExpression; 
      CompiledExpression<real_v(const real_t*, const complex_v*)>  m_probExpression; 
  };
} // namespace AmpGen

#endif
