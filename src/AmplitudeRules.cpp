#include "AmpGen/AmplitudeRules.h"

#include <memory.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <numeric>

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitExpression.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;
using namespace std::complex_literals; 

Coupling::Coupling(MinuitParameter* re, MinuitParameter* im) :
  m_re(re),
  m_im(im)
{
  auto tokens = split( re->name(), '_' );
  if ( tokens.size() == 3 ) {
    m_prefix = tokens[0];
    m_name   = tokens[1];
  } else if ( tokens.size() == 2 ) {
    m_name   = tokens[0];
  }
  else {
    ERROR("Ill-formed decay descriptor: " << m_name );
  }
  m_particle = Particle(m_name);
  coordinateType coord = NamedParameter<coordinateType>("CouplingConstant::Coordinates", coordinateType::cartesian);
  angType degOrRad     = NamedParameter<angType>("CouplingConstant::AngularUnits"      , angType::rad);
  m_isCartesian = true; 
  
  if( coord == coordinateType::polar ) m_isCartesian = false; 
  
  if ( coord == coordinateType::Invalid){
    FATAL("Coordinates for coupling constants must be either cartesian or polar");
  } 
  if ( degOrRad == angType::deg) m_sf = M_PI / 180; 
  
  if ( degOrRad == angType::Invalid ){
    FATAL("TotalCoupling::AngularUnits must be either rad or deg");
  } 
}

Coupling::Coupling(MinuitExpression* expression) : 
  m_name(expression->name()),
  m_expr(expression),
  m_particle(m_name){}

AmplitudeRules::AmplitudeRules( const MinuitParameterSet& mps )
{
  for ( auto& it_re : mps ) {
    auto& name = it_re->name(); 
    DEBUG("Attempting to parse: " << it_re->name() );
    if ( name.find("_Re") != std::string::npos ){
      auto it_im = mps.find(replaceAll( name, "_Re","_Im") );
      if( it_im == nullptr ){
        ERROR("Cannot find matching imaginary part / phase for: " <<  it_re->name() );
        continue; 
      }
      if( ! Particle::isValidDecayDescriptor( name.substr(0, name.find("_Re") ) ) ) continue; 
      Coupling p(it_re, it_im);
      m_rules[p.head()].emplace_back(p);
    }
    else if( name.find("_Im") == std::string::npos ){
      bool isCoupling = Particle::isValidDecayDescriptor( it_re->name() );
      if( isCoupling ){
        MinuitExpression* expression = dynamic_cast<MinuitExpression*>( it_re );
        DEBUG("Constructing: " << expression << " " << it_re->name() );
        if( expression != nullptr ){
          Coupling p(expression);
          m_rules[p.head()].emplace_back(p);
        } 
      } 
    }
  }
}

TotalCoupling::TotalCoupling(const TotalCoupling& other, const Coupling& pA) : 
  couplings(other.couplings)
{
  couplings.emplace_back(pA);
}


bool AmplitudeRules::hasDecay(const std::string& head) 
{ 
  return m_rules.find(head) != m_rules.end(); 
}

std::vector<Coupling> AmplitudeRules::rulesForDecay(const std::string& head, const std::string&  prefix)
{
  if(!hasDecay(head)) return std::vector<Coupling>();
  if( prefix == "" )return m_rules[head];
  std::vector<Coupling> rt = m_rules[head];
  rt.erase( std::remove_if( std::begin(rt), std::end(rt), [&prefix](auto& p){ return p.prefix() != prefix; } ), rt.end() );
  return rt;
}

const std::map<std::string, std::vector<Coupling>>& AmplitudeRules::rules() const
{ 
  return m_rules;
}

EventType Coupling::eventType() const
{
  Particle particle( m_name );
  std::vector<std::string> particleNames = { particle.name() };
  std::vector<std::shared_ptr<Particle>> fs = particle.getFinalStateParticles();
  std::stable_sort( fs.begin(), fs.end(), []( const auto& A, const auto& B ) { return *A < *B; } );
  std::transform( fs.begin(), fs.end(), std::back_inserter(particleNames), [](const auto& p) -> std::string { return p->name() ; } );
  return EventType( particleNames );
}

TotalCoupling::TotalCoupling(const Coupling& pA)
{
  couplings.emplace_back(pA);  
}

std::complex<double> Coupling::operator()() const 
{
  return m_expr != nullptr ? m_expr->getVal() : ( m_isCartesian ? complex_t( m_re->mean(), m_im->mean() ) : m_re->mean() * exp( 1i* m_sf * m_im->mean() ) ); 
}

Expression Coupling::to_expression() const 
{
  return m_expr != nullptr ? m_expr->expression() : ( m_isCartesian ? ComplexParameter(Parameter(m_re->name()), Parameter(m_im->name())) : Parameter( m_re->name() ) * fcn::exp( 1i * m_sf * Parameter(m_im->name()) ) );
}

std::complex<double> TotalCoupling::operator()() const
{
  return std::accumulate( couplings.begin(), couplings.end(), complex_t(1,0), [](const auto& prod, const auto& coupling ){ return prod * coupling() ; } );
}

Expression TotalCoupling::to_expression() const
{
  return std::accumulate( couplings.begin(), couplings.end(), Expression(1), [](const auto& prod, const auto& coupling ){ return prod * coupling.to_expression(); } );
}

void TotalCoupling::print() const
{
  INFO( this->operator()() << " -> "  );
  for( const auto& coupling :*this ) INFO( coupling.to_expression() << " = " << coupling() );
}

std::vector<std::pair<Particle, TotalCoupling>> AmplitudeRules::getMatchingRules(const EventType& type, const std::string& prefix )
{
  auto rules        = rulesForDecay( type.mother() );
  std::vector<std::pair<Particle, TotalCoupling>> rt; 
  for ( auto& rule : rules ) {
    if ( rule.prefix() != prefix ) continue;
    std::vector<std::pair<Particle, TotalCoupling>> tmpParticles;
    auto fs = type.finalStates();
    tmpParticles.emplace_back( Particle( rule.name(), fs ), TotalCoupling(rule) );
    do {
      std::vector<std::pair<Particle, TotalCoupling>> newTmpParticles;
      for ( auto& particleWithTotalCoupling : tmpParticles ) {
        auto protoParticle    = particleWithTotalCoupling.first;
        auto coupling         = particleWithTotalCoupling.second;
        auto protoFinalStates = protoParticle.getFinalStateParticles();
        if ( protoFinalStates.size() == type.size() ) {
          rt.emplace_back( particleWithTotalCoupling );
          continue;
        }
        std::string nameToExpand = protoParticle.uniqueString();
        for ( auto& ifs : protoFinalStates ) {
          auto expandedRules = rulesForDecay( ifs->name() ); /// get rules for decaying particle
          if ( expandedRules.size() == 0 ) continue;
          for ( auto& subTree : expandedRules ) {
            auto expanded_amplitude = replaceAll( nameToExpand, ifs->name(), subTree.name() );
            auto fs2                = type.finalStates();
            newTmpParticles.emplace_back( Particle( expanded_amplitude, fs2 ), TotalCoupling( coupling, subTree) );
          }
          break; // we should only break if there are rules to be expanded ...
        }
      }
      tmpParticles = newTmpParticles;
    } while ( tmpParticles.size() != 0 );
  }
  rt.erase( std::remove_if( std::begin(rt), std::end(rt), [](auto& p){ return !p.first.isStateGood(); } ), rt.end() );
  auto end = std::end(rt);
  for (auto it = rt.begin(); it != end; ++it) {
    auto dd = it->first.decayDescriptor();
    end = std::remove_if(it + 1, end, [dd](auto p){ return p.first.decayDescriptor() == dd;} );
  }
  rt.erase(end, rt.end());
  return rt;
}

bool TotalCoupling::isFixed() const
{
  return std::all_of(begin(), end(), [](auto& c){ return c.x()->isFixed() && c.y()->isFixed() ; } );
}

bool TotalCoupling::contains( const std::string& label ) const 
{
  return std::any_of(begin(), end(), [&label](auto& c){ return c.name().find(label) != std::string::npos ; } );
}

