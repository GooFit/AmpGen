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

#include "TMatrixD.h"

namespace AmpGen
{
  class LinearErrorPropagator;
  class MinuitParameter;
  class MinuitParameterSet;
  class FitFraction;
  class Event;
  class EventList;
  class MinuitProxy; 

  class PolarisedSum
  {
    public: 
      PolarisedSum() = default; 
      PolarisedSum( const EventType& eventType, AmpGen::MinuitParameterSet& mps, const std::string& prefix="" );
      void prepare();
      void setEvents( AmpGen::EventList& events );
      void setMC( AmpGen::EventList& events );
      void reset( const bool& flag = false );
      void debug(const AmpGen::Event& event );
      void debug_norm(); 
      void setWeight( MinuitParameter* param );
      double getWeight() const;
      void calculateNorms(const std::vector<bool>& hasChanged); 
      void generateSourceCode( const std::string& fname, const double& normalisation = 1, bool add_mt = false );
      void build_probunnormalised();
      Expression probExpression( const Tensor& T_matrix, const std::vector<Expression>& p, DebugSymbols* = nullptr ) const; 
      size_t size() const ;  
      real_t norm() const;
      complex_t norm(const size_t& i, const size_t& j, Integrator<18>* integ = nullptr ); 
      real_t operator()(const AmpGen::Event& event) const { return prob_unnormalised(event) ;}
      real_t prob_unnormalised( const AmpGen::Event& evt ) const;
      real_t prob(const AmpGen::Event& evt ) const ;
      real_t getValNoCache( const AmpGen::Event& evt ) ;
      std::vector<FitFraction> fitFractions( const LinearErrorPropagator& prop );
      std::vector<TransitionMatrix<std::vector<complex_t>>> matrixElements() const;
      void transferParameters(); 
      Tensor transitionMatrix();
      TransitionMatrix<std::vector<complex_t>> operator[](const size_t& i) const { return m_matrixElements[i] ; } 

    private: 
      size_t                        m_nCalls      = {0};
      real_t                        m_norm        = {1};
      EventList*                    m_events      = {nullptr};
      MinuitParameter*              m_weightParam = {nullptr}; 
      const MinuitParameterSet*     m_mps         = {nullptr};
      double                        m_weight      = {1}; 
      std::vector<MinuitProxy>      m_pVector     = {}; 
      bool                          m_verbosity   = {0};
      bool                          m_debug       = {0};
      Integrator<18>                m_integrator;
      std::vector<Bilinears>        m_norms;
      std::vector<std::vector<int>> m_polStates; 
      EventType                     m_eventType;
      std::string                   m_prefix      = "";
      std::vector<complex_t>        m_rho;
      std::vector<size_t>           m_integIndex; 
      std::vector<TransitionMatrix<std::vector<complex_t>>>        m_matrixElements;  
      CompiledExpression<real_t, const real_t*, const complex_t*> m_probExpression; 
      AmplitudeRules                m_rules;  
      std::vector<std::vector<int>> polarisationOuterProduct( const std::vector<std::vector<int>>& A, const std::vector<int>& B ) const;
      std::vector<int> polarisations( const std::string& name ) const ;
  };
} // namespace AmpGen

#endif
