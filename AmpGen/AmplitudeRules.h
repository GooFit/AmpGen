#ifndef AMPGEN_AMPLITUDERULES_H
#define AMPGEN_AMPLITUDERULES_H

#include <complex>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Particle.h"

namespace AmpGen
{
  class AmplitudeRules;
  class MinuitParameter;
  class MinuitParameterSet;
  struct Coupling;

  class AmplitudeRule
  {
  private:
    std::string m_prefix;
    std::string m_name;
    MinuitParameter* m_re;
    MinuitParameter* m_im;
    std::string m_head;
    bool m_isGood;
    std::vector<unsigned int> m_daughterProcesses;

  public:
    AmplitudeRule( const std::string& reName, std::map<std::string, MinuitParameter*>& mapping );
    std::string name() const { return m_name; }
    std::string prefix() const { return m_prefix; }
    std::string head() const { return m_head; }
    EventType eventType() const;
    friend struct Coupling;
    friend class AmplitudeRules;
    Coupling makeCoupling( bool isCartan = true );
  };

  struct Coupling {
    Coupling( const Coupling& other, const AmplitudeRule& pA, bool _isCartesian = true ) : couplings( other.couplings )
    {
      couplings.emplace_back( pA.m_re, pA.m_im );
      isCartesian = _isCartesian;
    }
    Coupling() { isCartesian = true; };
    std::vector<std::pair<MinuitParameter*, MinuitParameter*>> couplings;
    bool isCartesian;
    std::complex<double> operator()() const;
    Expression to_expression() const;
    void print() const;
    std::pair<MinuitParameter*, MinuitParameter*> operator[]( const size_t& index ) { return couplings[index]; }
  };

  class AmplitudeRules
  {
  private:
    std::map<std::string, std::vector<AmplitudeRule>> m_rules;

  public:
    AmplitudeRules() = default; 
    AmplitudeRules( MinuitParameterSet& mps );
    std::vector<AmplitudeRule> rulesForDecay( const std::string& head );
    bool hasDecay( const std::string& head );
    std::map<std::string, std::vector<AmplitudeRule>> rules();
    

    std::vector< std::pair< AmpGen::Particle, Coupling > > getMatchingRules( 
       const AmpGen::EventType& type, const std::string& prefix, const bool& useCartesian );
  };
} // namespace AmpGen

#endif
