#ifndef AMPGEN_AMPLITUDERULES_H
#define AMPGEN_AMPLITUDERULES_H

#include <stddef.h>
#include <complex>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Event.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ExpressionParser.h"
#if ENABLE_AVX 
  #include "AmpGen/EventListSIMD.h"
#endif

namespace AmpGen
{
  class MinuitParameter;
  class MinuitExpression;
  class MinuitParameterSet;

  class Coupling
  {
    public:
      Coupling(MinuitParameter* re, MinuitParameter* im);
      Coupling(MinuitExpression* expression);
      std::string name() const { return m_name; }
      std::string head() const { return m_particle.name(); }
      std::string prefix() const { return m_prefix; }
      EventType eventType() const;
      MinuitParameter* x() const { return m_re; }
      MinuitParameter* y() const { return m_im; }
      complex_t operator()() const; 
      Expression to_expression() const; 
      const Particle& particle() const { return m_particle ; }
      bool isCartesian() const { return m_isCartesian; }
    private:
      std::string       m_prefix = {""};
      std::string       m_name   = {""};
      MinuitParameter*  m_re     = {nullptr};
      MinuitParameter*  m_im     = {nullptr};
      MinuitExpression* m_expr   = {nullptr};
      Particle          m_particle;
      bool              m_isCartesian = {true};
      double            m_sf          = {1};
  };
  
  class TotalCoupling
  {
    public:
      TotalCoupling() = default; 
      TotalCoupling( const TotalCoupling& other, const Coupling& pA);
      TotalCoupling( const Coupling& pA);
      std::complex<double> operator()() const;
      Expression to_expression() const;
      void print() const;
      Coupling operator[]( const size_t& index ) { return couplings[index]; }
      bool isFixed() const; 
      bool contains( const std::string& name ) const;
      size_t size() const { return couplings.size(); }
      std::vector<Coupling>::const_iterator begin() const { return couplings.begin() ; }
      std::vector<Coupling>::const_iterator   end() const { return couplings.end() ; }
    private:
      std::vector<Coupling> couplings;
  };
  
  class AmplitudeRules
  {
    public:
      AmplitudeRules() = default; 
      AmplitudeRules( const MinuitParameterSet& mps );
      std::vector<Coupling> rulesForDecay(const std::string& head, const std::string& prefix="");
      bool hasDecay( const std::string& head );
      const std::map<std::string, std::vector<Coupling>>& rules() const;
      std::vector<std::pair<Particle, TotalCoupling>> getMatchingRules( 
          const EventType& type, const std::string& prefix="" );
      std::vector<Coupling> processesThatProduce(const Particle& particle) const; 

    private:
      std::map<std::string, std::vector<Coupling>> m_rules;
  };

  template <class RT> struct TransitionMatrix : public CompiledExpression<RT(const real_t*, const float_v*)>  
  {
    using amp_type = CompiledExpression<RT(const real_t*, const float_v*)>; 
    TransitionMatrix() = default;
    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const amp_type& amp) :
          amp_type(amp),
          decayTree(dt), 
          coupling(coupling) {}

    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const MinuitParameterSet& mps,
                     const std::map<std::string, unsigned>& evtFormat, 
                     const bool& debugThis=false) :
      amp_type(Particle(dt).getExpression(debugThis ? &db : nullptr ), dt.decayDescriptor(), evtFormat, db, &mps ),
      decayTree(dt),
      coupling(coupling) {}

    #if ENABLE_AVX
    const RT operator()(const Event& event) const { return amp_type::operator()(EventListSIMD::makeEvent(event).data()); }
    void debug( const Event& event ) const {               amp_type::debug(EventListSIMD::makeEvent(event).data() ) ; } 
    
    #else
    const RT operator()(const Event& event) const { return amp_type::operator()(event.address()) ; }
    void debug( const Event& event )        const {        amp_type::debug(event.address()) ; }
    #endif
    template <class... arg_types> auto operator()(arg_types... args ) const { return amp_type::operator()(args...) ; }
    
    const RT operator()(const float_v* t) const     { return amp_type::operator()(t) ; }
    void debug( const float_v* t )        const     {        amp_type::debug(t) ; } 
    const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  

    Particle                                            decayTree;
    TotalCoupling                                       coupling;
    complex_t                                           coefficient;
    DebugSymbols                                        db; 
    bool                                                workToDo    = {false};
  };
 
  template <class RT> std::vector<size_t> processIndex(const std::vector<TransitionMatrix<RT>>& tm, const std::string& label)
  {
    std::vector<size_t> indices;
    for ( size_t i = 0; i < tm.size(); ++i ) {
      if ( tm[i].coupling.contains(label) ) indices.push_back(i);
    }
    return indices;
  }
   
  template <class RT> size_t findIndex(const std::vector<TransitionMatrix<RT>>& tm, const std::string& decayDescriptor)
  {
    for ( size_t i = 0; i < tm.size(); ++i ) {
      if ( tm[i].decayDescriptor() == decayDescriptor ) return i;
    }
    ERROR( "Component " << decayDescriptor << " not found" );
    return 999;
  }

  template <class RT>   
  std::vector<size_t> findIndices(const std::vector<TransitionMatrix<RT>>& tm, const std::string& decayDescriptor)
  {
    std::vector<size_t> rt; 
    for ( size_t i = 0; i < tm.size(); ++i ) {
      if ( tm[i].decayDescriptor().find(decayDescriptor) != std::string::npos ) rt.push_back(i);
    }
    return rt;
  }
  
  template <> struct TransitionMatrix<void> : public CompiledExpression<void(complex_v*, const size_t&, const real_t*, const float_v*)>
  {
    using amp_type    = CompiledExpression<void(complex_v*, const size_t&, const real_t*, const float_v*)>;
    TransitionMatrix() = default;
    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const amp_type& amp) : 
          amp_type(amp),
          decayTree(dt), 
          coupling(coupling) {}

    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const MinuitParameterSet& mps,
                     const std::map<std::string, unsigned>& evtFormat, 
                     const bool& debugThis=false) :
      amp_type(Particle(dt).getExpression(debugThis ? &db : nullptr ), dt.decayDescriptor(), evtFormat, db, &mps ),
      decayTree(dt),
      coupling(coupling)
      { use_rto();}

    const std::vector<complex_v> operator()(const Event& event) const 
    { 
      std::vector<complex_v> rt(size); 
      #if ENABLE_AVX 
      amp_type::operator()(rt.data(), 1, externBuffer().data(), EventListSIMD::makeEvent(event).data());
      #else
      amp_type::operator()(rt.data(), 1, externBuffer().data(), event.address()); 
      #endif
      return rt;
    }
    template <class... arg_types> auto operator()(arg_types... args ) const { return amp_type::operator()(args...) ; }
    #if ENABLE_AVX
    void debug( const Event& event ) const { amp_type::debug(EventListSIMD::makeEvent(event).data() ) ; } 
    #else 
    void debug( const Event& event ) const { amp_type::debug(event.address()) ; }
    #endif 
    const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  

    Particle                                            decayTree;
    TotalCoupling                                       coupling;
    complex_t                                           coefficient;
    DebugSymbols                                        db; 
    bool                                                workToDo    = {false};
    unsigned                                            size        = {0};
  };
 
} // namespace AmpGen

#endif
