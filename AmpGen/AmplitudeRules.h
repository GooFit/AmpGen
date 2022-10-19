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
      Coupling(const Particle& particle, double f);
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
      static const AmplitudeRules* get(); 
      static AmplitudeRules* create( const MinuitParameterSet& mps, const std::string& prefix=""); 
      AmplitudeRules() = default;
      AmplitudeRules( const MinuitParameterSet& mps, const std::string& prefix="" );
      std::vector<Coupling> rulesForDecay(const std::string& head, const std::string& prefix="") const;
      bool hasDecay( const std::string& head ) const;
      const std::map<std::string, std::vector<Coupling>>& rules() const;
      std::vector<std::pair<Particle, TotalCoupling>> getMatchingRules( 
          const EventType& type, const std::string& prefix="" );
      std::vector<Coupling> processesThatProduce(const Particle& particle) const; 

      std::vector<std::pair<Particle, TotalCoupling>> expand( const Coupling& coupling ) const; 
      void add_rule( const Particle& p , double coupling ) ;
    private:
      std::map<std::string, std::vector<Coupling>> m_rules;
      static AmplitudeRules* gAmplitudeRules;
  };
 
} // namespace AmpGen

#endif
