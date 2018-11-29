#include "AmpGen/FastIncoherentSum.h"

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

namespace AmpGen
{
  class EventType;
  class MinuitParameterSet;
} // namespace AmpGen

using namespace AmpGen;

FastIncoherentSum::FastIncoherentSum( const EventType& finalStates, AmpGen::MinuitParameterSet& mps,
                                      const std::string& prefix )
    : FastCoherentSum( finalStates, mps, prefix )
{

  m_normalisations.resize( size(), 1 );
}

double FastIncoherentSum::norm() const
{
  double norm( 0 ); // (0,0);
  for ( unsigned int i = 0; i < size(); ++i ) {
    double val = m_normalisations.get( i, 0 ).real() * std::norm( m_matrixElements[i].coefficient );
    norm += val;
  }
  return norm; //.real();
}

double FastIncoherentSum::norm( const Bilinears& norms ) const
{
  double norm( 0 ); // (0,0);
  for ( unsigned int i = 0; i < size(); ++i ) {
    double val = norms.get( i, 0 ).real() * std::norm( m_matrixElements[i].coefficient );
    norm += val;
  }
  return norm; //.real();
}

void FastIncoherentSum::prepare()
{
  if ( m_weightParam != nullptr ) m_weight = m_weightParam->mean();

  if ( m_isConstant && m_prepareCalls != 0 ) return;

  transferParameters(); /// move everything to the "immediate" cache ///
  bool isReady = true; 
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) {
    auto& pdf = m_matrixElements[i].pdf;
    if( m_prepareCalls == 0 ) isReady &= pdf.isReady();
    pdf.prepare();
    if ( m_prepareCalls != 0 && !pdf.hasExternalsChanged() ) continue;

    if ( m_prepareCalls == 0 && m_events != nullptr )
      m_matrixElements[i].addressData                                             = m_events->registerExpression( pdf );

    if ( m_events != nullptr ) m_events->updateCache( pdf, m_matrixElements[i].addressData );
    
    if ( m_prepareCalls == 0 && m_integrator.isReady() ) m_integrator.prepareExpression( pdf );

    pdf.resetExternals();
  }
  m_integrator.flush();

  m_prepareCalls++;
  m_norm = norm(); /// update normalisation
}

std::vector<FitFraction> FastIncoherentSum::fitFractions( const LinearErrorPropagator& linProp )
{
  std::vector<FitFraction> outputFractions;
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) {

    IFFCalculator calc;
    calc.fcs   = this;
    calc.index = i;

    outputFractions.emplace_back( m_matrixElements[i].decayTree->uniqueString(), calc(), linProp.getError( calc ) );
  }

  for ( auto& p : outputFractions ) {
    INFO( std::setw( 100 ) << p.name() << " " << std::setw( 7 ) << p.val() << " ± " << p.err() );
  }
  return outputFractions;
}
