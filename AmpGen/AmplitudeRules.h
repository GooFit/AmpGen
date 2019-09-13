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

namespace AmpGen
{
  class AmplitudeRules;
  class MinuitParameter;
  class MinuitParameterSet;
  class CouplingConstant;

  class AmplitudeRule
  {
    public:
      AmplitudeRule(MinuitParameter* re, MinuitParameter* im);
      std::string name() const { return m_name; }
      std::string head() const { return m_particle.name(); }
      std::string prefix() const { return m_prefix; }
      EventType eventType() const;
      friend class CouplingConstant;
      friend class AmplitudeRules;

    private:
      std::string      m_prefix = {""};
      std::string      m_name   = {""};
      MinuitParameter* m_re     = {nullptr};
      MinuitParameter* m_im     = {nullptr};
      Particle         m_particle;
  };
  
  class AmplitudeRules
  {
    public:
      AmplitudeRules() = default; 
      AmplitudeRules( const MinuitParameterSet& mps );
      std::vector<AmplitudeRule> rulesForDecay(const std::string& head, const std::string& prefix="");
      bool hasDecay( const std::string& head );
      std::map<std::string, std::vector<AmplitudeRule>> rules();
      std::vector<std::pair<Particle, CouplingConstant>> getMatchingRules( 
          const EventType& type, const std::string& prefix="" );
      std::vector<AmplitudeRule> processesThatProduce(const Particle& particle) const; 

    private:
      std::map<std::string, std::vector<AmplitudeRule>> m_rules;
  };

  class CouplingConstant 
  {
    public:
      CouplingConstant() = default; 
      CouplingConstant( const CouplingConstant& other, const AmplitudeRule& pA);
      CouplingConstant( const AmplitudeRule& pA);
      std::complex<double> operator()() const;
      Expression to_expression() const;
      void print() const;
      std::pair<MinuitParameter*, MinuitParameter*> operator[]( const size_t& index ) { return couplings[index]; }
      bool isFixed() const; 
      bool contains( const std::string& name ) const;
      std::vector<std::pair<MinuitParameter*, MinuitParameter*>> couplings;

    private:
      bool isCartesian = {true};
      double sf        = {1};
  };

  template <class RT> struct TransitionMatrix 
  {
    TransitionMatrix() = default;
    TransitionMatrix(const Particle& dt, 
                     const CouplingConstant& coupling, 
                     const CompiledExpression<RT, const real_t*, const real_t*> & amp) : 
          decayTree(dt), 
          coupling(coupling), 
          amp(amp) {}

    TransitionMatrix(Particle& dt, 
                     const CouplingConstant& coupling, 
                     const MinuitParameterSet& mps,
                     const std::map<std::string,size_t>& evtFormat, 
                     const bool& debugThis=false) :
      decayTree(dt),
      coupling(coupling),
      amp(decayTree.getExpression(debugThis ? &db : nullptr ), decayTree.decayDescriptor(), evtFormat, db, &mps ) {}

    const RT operator()(const Event& event) const { return amp(event.address() ); }
    const RT operator()(const Event& event, const size_t& cacheOffset) const { return amp(event.address() + cacheOffset); }
    const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  

    Particle                                            decayTree;
    CouplingConstant                                    coupling;
    complex_t                                           coefficient;
    DebugSymbols                                        db; 
    CompiledExpression<RT,const real_t*,const real_t*>  amp; 
    size_t                                              addressData = {999};
  };
 
  template <class RT>  
  std::vector<size_t> processIndex(const std::vector<TransitionMatrix<RT>>& tm, const std::string& label)
  {
    std::vector<size_t> indices;
    for ( size_t i = 0; i < tm.size(); ++i ) {
      if ( tm[i].coupling.contains(label) ) indices.push_back(i);
    }
    return indices;
  }
   
  template <class RT>
  size_t findIndex(const std::vector<TransitionMatrix<RT>>& tm, const std::string& decayDescriptor)
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

} // namespace AmpGen

#endif
