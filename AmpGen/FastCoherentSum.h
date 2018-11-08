#ifndef AMPGEN_FASTCOHERENTSUM_H
#define AMPGEN_FASTCOHERENTSUM_H

#include <memory.h>
#include <stddef.h>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Types.h"
#include "AmpGen/Event.h"

namespace AmpGen
{
  class LinearErrorPropagator;
  class MinuitParameter;
  class MinuitParameterSet;
  class FitFraction;
  class Particle;

  template <class RT>
  struct TransitionMatrix {
    std::shared_ptr<Particle>                           decayTree;
    Coupling                                            coupling;
    complex_t                                           coefficient;
    CompiledExpression<RT,const real_t*,const real_t*>  pdf; 
    size_t                                              addressData;
    size_t                                              addressInt;
    const RT operator()( const Event& event ) const { return pdf(event.address() ); }
    TransitionMatrix(){};
    TransitionMatrix( const std::shared_ptr<Particle>& dt, 
                      const Coupling& coup, 
                      const CompiledExpression<RT, const real_t*, const real_t*> & _pdf ) : 
          decayTree( dt ), 
          coupling( coup ), 
          pdf( _pdf ), 
          addressData( 999 ), 
          addressInt( 999 )
    {
    }
    std::vector<MinuitParameter*> getDependencies() const;
  };

  class FastCoherentSum
  {
  protected:
    std::vector<TransitionMatrix<complex_t>> m_matrixElements; ///< Vector of (expanded) matrix elements
    Bilinears m_normalisations;                  ///< Normalisation integrals
    AmplitudeRules m_protoAmplitudes;            ///< Proto amplitudes from user rule-set
    Integrator<10> m_integralDispatch;           ///< Integral dispatch tool (with default unroll = 10) 
    TransitionMatrix<complex_t> m_total;         ///< Total Matrix Element 
    EventList* m_events = {nullptr};             ///< Data events to evaluate PDF on
    EventList* m_sim = {nullptr};                ///< Integration events to evalute PDF on
    EventType m_evtType;                         ///< Final state for this amplitude
    MinuitParameter* m_weightParam = {nullptr};  ///< Weight parameter (i.e. the normalised yield)
    int m_prepareCalls             = {0};        ///< Number of times prepare has been called
    int m_lastPrint                = {0};        ///< Last time verbose PDF info was printed
    int m_printFreq                = {0};        ///< Frequency to print verbose PDF info
    double m_weight                = {1};        ///< Weight number (i.e. the normalised yield)
    double m_norm                  = {0};        ///< Normalisation integral
    std::string m_prefix           = {""};       ///< Prefix for matrix elements
    bool m_stateIsGood             = {true};     ///< Flag for the state being good
    bool m_isConstant              = {false};    ///< Flag for a constant PDF
    bool m_dbThis                  = {false};    ///< Flag to generate amplitude level debugging
    bool m_verbosity               = {false};    ///< Flag for verbose printing

    void addMatrixElement( std::pair<Particle, Coupling>& particleWithCoupling, const MinuitParameterSet& mps );
    bool isFixedPDF(const MinuitParameterSet& mps) const;

  public:
    FastCoherentSum( const EventType& type, AmpGen::MinuitParameterSet& mps, const std::string& prefix = "" );

    real_t operator()( const Event& evt ) const { return prob( evt ); }

    AmplitudeRules protoAmplitudes() { return m_protoAmplitudes; }
    std::vector<TransitionMatrix<complex_t>> matrixElements() { return m_matrixElements; }

    std::vector<unsigned int> processIndex( const std::string& label ) const;
    std::string getParentProcess( const std::string& label ) const;

    unsigned int getPdfIndex( const std::string& name ) const;
    void PConjugate();

    TransitionMatrix<complex_t> operator[]( const unsigned int& index ) { return m_matrixElements[index]; }
    const TransitionMatrix<complex_t> operator[]( const unsigned int& index ) const { return m_matrixElements[index]; }

    std::string prefix() const { return m_prefix; }

    unsigned int size() const { return m_matrixElements.size(); }
    double getWeight() const { return m_weight; }
    void setWeight( const double& weight ) { m_weight = weight; }
    void setWeight( MinuitParameter* param ) { m_weightParam = param; }

    void makeTotalExpression();
    void reset( bool resetEvents = false );
    void setEvents( EventList& list );
    void setMC( EventList& sim );

    real_t norm( const Bilinears& norms ) const;
    real_t norm() const;
    bool isStateGood() const { return m_stateIsGood; }
    complex_t norm( const unsigned int& x, const unsigned int& y ) const;
    void transferParameters();
    void preprepare();
    void prepare();
    void printVal( const Event& evt, bool isSim = false );
    void updateNorms( const std::vector<unsigned int>& changedPdfIndices );
    
    std::vector<unsigned int> cacheAddresses( const EventList& evts ) const
    {
      std::vector<unsigned int> addresses;
      for ( auto& mE : m_matrixElements ) {
        addresses.push_back( evts.getCacheIndex( mE.pdf ) );
      }
      return addresses;
    }

    complex_t getVal( const Event& evt ) const
    {
      complex_t value( 0., 0. );
      for ( auto& mE : m_matrixElements ) {
        value += mE.coefficient * evt.getCache( mE.addressData );
      }
      return value;
    }

    complex_t getVal( const Event& evt, const std::vector<unsigned int>& cacheAddresses ) const
    {
      complex_t value( 0., 0. );
      for ( unsigned int i = 0; i < m_matrixElements.size(); ++i )
        value += m_matrixElements[i].coefficient * evt.getCache( cacheAddresses[i] );
      return value;
    }
    std::complex<double> getValNoCache( const Event& evt ) const;

    double prob( const Event& evt ) const { return m_weight * std::norm( getVal( evt ) ) / m_norm; }
    double prob_unnormalised( const Event& evt ) const { return std::norm( getVal( evt ) ); }

    void debug( const Event& evt, const std::string& nameMustContain="");
    std::vector<FitFraction> fitFractions( const LinearErrorPropagator& linProp );

    void generateSourceCode( const std::string& fname, const double& normalisation = 1, bool add_mt = false );

    double getNorm( const Bilinears& normalisations );
    
    double get_norm() const { return m_norm ; }
    std::map<std::string, std::vector<unsigned int>> getGroupedAmplitudes();
    void resync();
    Bilinears norms() const { return m_normalisations ; }
  };
} // namespace AmpGen

#endif
