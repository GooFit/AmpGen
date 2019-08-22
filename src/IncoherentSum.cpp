#include "AmpGen/IncoherentSum.h"

#include <memory.h>
#include <iomanip>
#include <memory>
#include <ostream>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/EventList.h"

using namespace AmpGen;

IncoherentSum::IncoherentSum( const EventType& finalStates, 
                              const MinuitParameterSet& mps,
                              const std::string& prefix ) : 
  CoherentSum( finalStates, mps, prefix )
{
  m_normalisations.resize( size(), 1 );
}

double IncoherentSum::norm() const
{ 
  return norm( m_normalisations ); 
}

double IncoherentSum::norm( const Bilinears& norms ) const
{
  double norm( 0 ); // (0,0);
  for ( unsigned int i = 0; i < size(); ++i ) {
    double val = norms.get( i, 0 ).real() * std::norm( m_matrixElements[i].coefficient );
    norm += val;
  }
  return norm; //.real();
}

void IncoherentSum::prepare()
{
  if ( m_weightParam != nullptr ) m_weight = m_weightParam->mean();
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  transferParameters();
  for ( auto& mE : m_matrixElements ) {
    auto& amp = mE.amp;
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    if ( m_prepareCalls == 0 && m_events != nullptr )
      mE.addressData = m_events->registerExpression( amp );
    if ( m_events != nullptr ) m_events->updateCache( amp, mE.addressData ); 
    if ( m_prepareCalls == 0 && m_integrator.isReady() ){
      m_integrator.prepareExpression( amp );
    }
    INFO( mE.addressData << " " << m_events->at(0).getCache(mE.addressData) );
    amp.resetExternals();
  }
  if( m_prepareCalls == 0 ){
    for( size_t i = 0 ; i < m_matrixElements.size(); ++i ){
      auto index = m_integrator.events().getCacheIndex( m_matrixElements[i].amp );
      m_integrator.queueIntegral( index, index, i, 0, &m_normalisations, false);
    }
    m_integrator.flush();
  }
  m_prepareCalls++;
  m_norm = norm();
  INFO( "norm = " << m_norm << " weight = " << m_weight );
}

std::vector<FitFraction> IncoherentSum::fitFractions( const LinearErrorPropagator& linProp )
{
  std::vector<FitFraction> outputFractions;
  std::vector<size_t> normSet(m_matrixElements.size());
  std::iota( std::begin(normSet), std::end(normSet), 0 );
  FitFractionCalculator<IncoherentSum> calc(this, normSet ); 
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) 
    calc.emplace_back(m_matrixElements[i].decayDescriptor(), std::vector<size_t>({i}));
  for ( auto& p : outputFractions ) INFO(p);
  return outputFractions;
}

double IncoherentSum::prob( const Event& evt ) const 
{ 
  return m_weight * prob_unnormalised(evt) / m_norm; 
}

double IncoherentSum::prob_unnormalised( const Event& evt ) const 
{ 
  double value( 0. );
  for ( auto& mE : m_matrixElements ) {
    value += std::norm( mE.coefficient * evt.getCache( mE.addressData ) );
  }
  return value;
}
