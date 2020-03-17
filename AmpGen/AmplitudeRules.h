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

  template <class RT> struct TransitionMatrix 
  {
    TransitionMatrix() = default;
    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const CompiledExpression<RT, const real_t*, const real_t*> & amp) : 
          decayTree(dt), 
          coupling(coupling), 
          amp(amp) {}

    TransitionMatrix(Particle& dt, 
                     const TotalCoupling& coupling, 
                     const MinuitParameterSet& mps,
                     const std::map<std::string, unsigned>& evtFormat, 
                     const bool& debugThis=false) :
      decayTree(dt),
      coupling(coupling),
      amp(decayTree.getExpression(debugThis ? &db : nullptr ), decayTree.decayDescriptor(), evtFormat, db, &mps ) {}

    const RT operator()(const Event& event) const { return amp(event.address() ); }
    const RT operator()(const Event& event, const size_t& cacheOffset) const { return amp(event.address() + cacheOffset); }
    const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  

    Particle                                            decayTree;
    TotalCoupling                                       coupling;
    complex_t                                           coefficient;
    DebugSymbols                                        db; 
    CompiledExpression<RT,const real_t*,const real_t*>  amp; 
    size_t                                              addressData = {999};
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
  
  template <> struct TransitionMatrix<void>
  {
    TransitionMatrix() = default;
    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const CompiledExpression<void, complex_t*, const real_t*, const real_t*> & amp) : 
          decayTree(dt), 
          coupling(coupling), 
          amp(amp) {}

    TransitionMatrix(const Particle& dt, 
                     const TotalCoupling& coupling, 
                     const MinuitParameterSet& mps,
                     const std::map<std::string, unsigned>& evtFormat, 
                     const bool& debugThis=false) :
      decayTree(dt),
      coupling(coupling),
      amp(decayTree.getExpression(debugThis ? &db : nullptr ), decayTree.decayDescriptor(), evtFormat, db, &mps ) { amp.use_rto();}

    const std::vector<complex_t> operator()(const Event& event) const { 
      std::vector<complex_t> rt(4); 
      amp(rt.data(), amp.externBuffer().data(), event.address() ); 
      return rt;
    }
    const std::vector<complex_t> operator()(const Event& event, const size_t& cacheOffset) const { 
      std::vector<complex_t> rt(4); 
      amp(rt.data(), amp.externBuffer().data(), event.address() + cacheOffset); 
      return rt;
    }
    const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  

    Particle                                            decayTree;
    TotalCoupling                                       coupling;
    complex_t                                           coefficient;
    DebugSymbols                                        db; 
    CompiledExpression<void, complex_t*, const real_t*, const real_t*>  amp; 
    size_t                                              addressData = {999};
    bool                                                workToDo    = {false};
  };
 
} // namespace AmpGen

#endif
