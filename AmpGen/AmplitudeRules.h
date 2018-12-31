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

namespace AmpGen
{
  class AmplitudeRules;
  class MinuitParameter;
  class MinuitParameterSet;
  class CouplingConstant;
  class Particle;

  class AmplitudeRule
  {
    public:
      AmplitudeRule( const std::string& reName, const std::map<std::string, MinuitParameter*>& mapping );
      std::string name() const { return m_name; }
      std::string prefix() const { return m_prefix; }
      std::string head() const { return m_head; }
      EventType eventType() const;
      friend class CouplingConstant;
      friend class AmplitudeRules;

    private:
      std::string m_prefix;
      std::string m_name;
      MinuitParameter* m_re;
      MinuitParameter* m_im;
      std::string m_head;
      bool m_isGood;
  };
  
  class AmplitudeRules
  {
    public:
      AmplitudeRules() = default; 
      AmplitudeRules( const MinuitParameterSet& mps );
      std::vector<AmplitudeRule> rulesForDecay( const std::string& head );
      bool hasDecay( const std::string& head );
      std::map<std::string, std::vector<AmplitudeRule>> rules();
      std::vector< std::pair< Particle, CouplingConstant > > getMatchingRules( 
          const EventType& type, const std::string& prefix, const bool& useCartesian );

    private:
      std::map<std::string, std::vector<AmplitudeRule>> m_rules;
  };

  class CouplingConstant 
  {
    public:
      CouplingConstant() = default; 
      CouplingConstant( const CouplingConstant& other, const AmplitudeRule& pA, bool isCartesian = true );
      CouplingConstant( const AmplitudeRule& pA, bool isCartesian = true );
      std::complex<double> operator()() const;
      Expression to_expression() const;
      void print() const;
      std::pair<MinuitParameter*, MinuitParameter*> operator[]( const size_t& index ) { return couplings[index]; }
      bool isFixed() const; 
      bool contains( const std::string& name ) const;
      void changeSign();
      std::vector<std::pair<MinuitParameter*, MinuitParameter*>> couplings;

    private:
      bool isCartesian = {true};
  };

  template <class RT> struct TransitionMatrix 
  {
    TransitionMatrix() = default;
    TransitionMatrix( const std::shared_ptr<Particle>& dt, 
                      const CouplingConstant& coup, 
                      const CompiledExpression<RT, const real_t*, const real_t*> & _pdf ) : 
          decayTree( dt ), 
          coupling( coup ), 
          pdf( _pdf ) {}
    const RT operator()( const Event& event ) const { return pdf(event.address() ); }
    
    std::shared_ptr<Particle>                           decayTree;
    CouplingConstant                                    coupling;
    complex_t                                           coefficient;
    CompiledExpression<RT,const real_t*,const real_t*>  pdf; 
    size_t                                              addressData = {999};
  };

} // namespace AmpGen

#endif
