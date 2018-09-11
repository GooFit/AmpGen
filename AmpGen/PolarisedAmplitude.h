#ifndef AMPGEN_POLARISEDAMPLITUDE_H
#define AMPGEN_POLARISEDAMPLITUDE_H

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Types.h"
#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Particle.h"
#include "AmpGen/FastCoherentSum.h"
#include "TMatrixD.h"

namespace AmpGen
{

  class CompiledExpressionBase;
  class LinearErrorPropagator;
  class MinuitParameter;
  class MinuitParameterSet;
  class FitFraction;
  class Particle;

  class PolarisedAmplitude
  {
    private: 
      std::vector< TransitionMatrix<std::vector<complex_t>>> m_matrixElements;  
      size_t                        m_size; 
      size_t                        m_nCalls;
      real_t                        m_norm;
      EventList*                    m_events;
      std::vector<MinuitProxy>      m_productionParameters; 
      Integrator<10>                m_integralDispatch;
      std::array<Bilinears,8>       m_norms;
      MinuitParameter*              m_weightParam; 
      MinuitParameterSet*           m_mps;
      double                        m_weight; 
      std::vector<std::vector<int>> m_polStates; 
      EventType                     m_eventType;
      CompiledExpression< double, const real_t*, const complex_t*> m_probExpression; 
    public: 
      PolarisedAmplitude( const EventType& eventType, AmpGen::MinuitParameterSet& mps, const std::string& prefix="" );
      void prepare();
      void setEvents( AmpGen::EventList& events );
      void setMC( AmpGen::EventList& events );
      void reset( const bool& flag = false );
      void debug(const AmpGen::Event& event ); 
      void setWeight( MinuitParameter* param ){ m_weightParam = param ; } 
      void calculateNorms(); 
      void generateSourceCode( const std::string& fname, const double& normalisation = 1, bool add_mt = false );
      void build_probunnormalised();
      Expression probExpression( const Tensor& T_matrix, const std::vector<Expression>& p ) const; 
      size_t size() const ;  
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
