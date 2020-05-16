#ifndef AMPGEN_COHERENTSUM_H
#define AMPGEN_COHERENTSUM_H

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
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Types.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Store.h"
#include "AmpGen/LiteSpan.h"
namespace AmpGen
{
  class LinearErrorPropagator;
  class MinuitParameterSet;
  class FitFraction;
  class Particle;

  /** @class CoherentSum
      @brief A coherent sum of amplitudes 

      An coherent sum of resonant contributions, where at a position in phase-space @f$\psi@f$, 
      the (unnormalised) probability density  is given by
      @f[ 
        \mathcal{P}(\psi) = \left| \sum_{i} g_i \mathcal{A}_i(\psi) \right|^2,    
      @f]
      where @f$\mathcal{P}(\psi)@f$ is the probability, @f$g_i@f$ is the coupling to an isobar channel, 
      and @f$\mathcal{A}_i(\psi)@f$ is the amplitude of the ith channel.
  */
  class CoherentSum
  {
  public:
    #if ENABLE_AVX
      using EventList_type = EventListSIMD; 
    #else 
      using EventList_type  = EventList; 
    #endif
    CoherentSum();
    CoherentSum( const EventType& type, const AmpGen::MinuitParameterSet& mps, const std::string& prefix = "" );
    virtual ~CoherentSum() = default; 

    AmplitudeRules protoAmplitudes() { return m_rules; }
    std::string prefix() const { return m_prefix; }
    
    auto operator[]( const size_t& index ) { return m_matrixElements[index]; }
    auto operator[]( const size_t& index ) const { return m_matrixElements[index]; }
    size_t size()                          const { return m_matrixElements.size(); }
    real_t getWeight()                     const { return m_weight; }
    real_t norm( const Bilinears& norms )  const;
    real_t norm() const;
    real_t getNorm( const Bilinears& normalisations );

    complex_t norm( const size_t& x, const size_t& y ) const;
    complex_t getVal( const Event& evt ) const;
    complex_t getValNoCache( const Event& evt ) const; 

    void transferParameters();
    void prepare();
    void printVal( const Event& evt );
    void updateNorms();
    void setWeight( MinuitProxy param ) { m_weight = param; }
    void makeTotalExpression();
    void reset( bool resetEvents = false );
    void setEvents( const EventList_type& list );
    void setMC( const EventList_type& sim );
    #if ENABLE_AVX
      void setEvents( const EventList& list) { 
        WARNING("Setting events from a AoS container, will need to make a copy");
        m_ownEvents = true; setEvents( *(new EventListSIMD(list)) ) ; } 
      void setMC    ( const EventList& list) { 
        WARNING("Setting integration events from a AoS container, will need to make a copy");
        setMC( *(new EventListSIMD(list)) ) ; } 
      double operator()(const double*, const unsigned) const; 
    #endif

    float_v operator()(const float_v*, const unsigned) const; 
    real_t  operator()(const Event& evt )              const { return m_weight*std::norm(getVal(evt))/m_norm; }

    void debug( const Event& evt, const std::string& nameMustContain="");
    void generateSourceCode( const std::string& fname, const double& normalisation = 1, bool add_mt = false );

    std::vector<FitFraction> fitFractions( const LinearErrorPropagator& linProp );
    auto matrixElements() const { return m_matrixElements; }

    std::map<std::string, std::vector<unsigned int>> getGroupedAmplitudes();
    Bilinears norms() const { return m_normalisations ; }
      
    std::function<real_t(const Event&)> evaluator(const EventList_type* = nullptr) const; 
    KeyedFunctors<double, Event> componentEvaluator(const EventList_type* = nullptr) const; 

  protected:
    std::vector<TransitionMatrix<complex_v>> m_matrixElements; ///< Vector of matrix elements
    Bilinears        m_normalisations;                         ///< Normalisation integrals
    AmplitudeRules   m_rules;                                  ///< Ruleset for the selected transition.
     
    Integrator       m_integrator;                             ///< Tool to calculate integrals 
    const EventList_type*  m_events       = {nullptr};               ///< Data events to evaluate PDF on
    Store<complex_v, Alignment::AoS> m_cache;                  ///< Store of intermediate values for the PDF calculation 

    bool             m_ownEvents    = {false};                 ///< Flag as to whether events are owned by this PDF or not  
    EventType        m_evtType;                                ///< Final state for this amplitude
    size_t           m_prepareCalls = {0};                     ///< Number of times prepare has been called
    size_t           m_lastPrint    = {0};                     ///< Last time verbose PDF info was printed
    size_t           m_printFreq    = {0};                     ///< Frequency to print verbose PDF info
    MinuitProxy      m_weight       = {nullptr, 1};            ///< Weight (i.e. the normalised yield)
    double           m_norm         = {1};                     ///< Normalisation integral
    bool             m_isConstant   = {false};                 ///< Flag for a constant PDF
    bool             m_dbThis       = {false};                 ///< Flag to generate amplitude level debugging
    bool             m_verbosity    = {false};                 ///< Flag for verbose printing
    std::string      m_objCache     = {""};                    ///< Directory that contains (cached) amplitude objects
    std::string      m_prefix       = {""};                    ///< Prefix for matrix elements
    const MinuitParameterSet* m_mps       = {nullptr};

    void addMatrixElement( std::pair<Particle, TotalCoupling>& particleWithCoupling, const MinuitParameterSet& mps );

  };
} // namespace AmpGen

#endif
