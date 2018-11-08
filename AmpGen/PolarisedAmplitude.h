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
#include "AmpGen/FastCoherentSum.h"
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

  class PolarisedAmplitude
  {
    private: 
      std::vector< TransitionMatrix<std::vector<complex_t>>> m_matrixElements;  
      size_t                        m_size = {0}; 
      size_t                        m_nCalls = {0};
      real_t                        m_norm = {1};
      EventList*                    m_events = { nullptr };
      std::vector<MinuitProxy>      m_productionParameters; 
      Integrator<10>                m_integralDispatch;
      std::array<Bilinears,8>       m_norms;
      MinuitParameter*              m_weightParam = {nullptr}; 
      MinuitParameterSet*           m_mps;
      double                        m_weight = {1}; 
      std::vector<std::vector<int>> m_polStates; 
      EventType                     m_eventType;
      CompiledExpression< real_t, const real_t*, const complex_t*> m_probExpression; 
     
      std::vector< std::vector< int > > polarisationOuterProduct( const std::vector< std::vector< int > >& A, const std::vector< int >& B ) const;
      std::vector<int> polarisations( const std::string& name ) const ;
    public: 
      PolarisedAmplitude( const EventType& eventType, AmpGen::MinuitParameterSet& mps, const std::string& prefix="" );
      void prepare();
      void setEvents( AmpGen::EventList& events );
      void setMC( AmpGen::EventList& events );
      void reset( const bool& flag = false );
      void debug(const AmpGen::Event& event );
      void debug_norm(); 
      void setWeight( MinuitParameter* param ){ m_weightParam = param ; } 
      void calculateNorms( const std::vector<size_t>& changedPdfIndices); 
      void generateSourceCode( const std::string& fname, const double& normalisation = 1, bool add_mt = false );
      void build_probunnormalised();
      Expression probExpression( const Tensor& T_matrix, const std::vector<Expression>& p ) const; 
      size_t size() const ;  
      complex_t TE( const AmpGen::Event& event, const size_t& x, const size_t& y );
      real_t norm() const;
      real_t prob_unnormalised( const AmpGen::Event& evt ) const;
      real_t prob(const AmpGen::Event& evt ) const ;
      real_t getValNoCache( const AmpGen::Event& evt ) ;
      std::vector<FitFraction> fitFractions( const LinearErrorPropagator& prop );
      std::vector<TransitionMatrix<std::vector<complex_t>>> matrixElements() const;
      void transferParameters(); 
      Tensor transitionMatrix();
  };
} // namespace AmpGen

#endif
