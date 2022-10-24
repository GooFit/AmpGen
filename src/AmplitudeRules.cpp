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

AmplitudeRules::AmplitudeRules(const MinuitParameterSet& mps)
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
      if( ! Particle::isValidDecayDescriptor(name.substr(0, name.find("_Re"))) ) continue;       
      Coupling p(it_re, it_im);
      m_rules[p.head()].emplace_back(p);
    }
    else if( name.find("_Im") == std::string::npos ){
      bool isCoupling = Particle::isValidDecayDescriptor( it_re->name() );
      if( !isCoupling ) continue; 
      MinuitExpression* expression = dynamic_cast<MinuitExpression*>( it_re );
      DEBUG("Constructing: " << expression << " " << it_re->name() );
      if( expression != nullptr ){
        Coupling p(expression);
        m_rules[p.head()].emplace_back(p);
      } 
    } 
  }
}

TotalCoupling::TotalCoupling(const TotalCoupling& other, const Coupling& pA) : 
  couplings(other.couplings)
{
  couplings.emplace_back(pA);
}


bool AmplitudeRules::hasDecay(const std::string& head) const 
{ 
  return m_rules.find(head) != m_rules.end(); 
}

std::vector<Coupling> AmplitudeRules::rulesForDecay(const std::string& head, const std::string&  prefix) const
{
  if(!hasDecay(head)) return std::vector<Coupling>();
  if( prefix == "" ) return m_rules.find(head)->second;
  std::vector<Coupling> rt = m_rules.find(head)->second;
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
  auto rules        = rulesForDecay( type.mother(), prefix );
  std::vector<std::pair<Particle, TotalCoupling>> rt; 
  for ( const auto& rule : rules ) {
    if ( rule.prefix() != prefix ) continue;
    auto expanded = expand( rule );
    for( const auto& [particle,coupling] : expanded )
    {
      auto p = Particle( particle.decayDescriptor(), type.finalStates() );
      if( p.isStateGood() && 
          std::find_if( rt.begin(), rt.end(), [p](const auto& x){ return x.first.decayDescriptor() == p.decayDescriptor(); } )
          == rt.end() 
        ) rt.emplace_back( p, coupling );
    }
  }
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

AmplitudeRules* AmplitudeRules::gAmplitudeRules = nullptr; 

AmplitudeRules* AmplitudeRules::create( const MinuitParameterSet& mps)
{
  if( gAmplitudeRules != nullptr )
  {
    WARNING("Recreating ruleset");
    delete gAmplitudeRules; 
  }
  gAmplitudeRules = new AmplitudeRules(mps);
  return gAmplitudeRules; 
}

const AmplitudeRules* AmplitudeRules::get() 
{
  if( gAmplitudeRules == nullptr ){ FATAL("No ruleset created"); } 
  return gAmplitudeRules; 
}

std::vector<std::pair<Particle, TotalCoupling>> AmplitudeRules::expand( const Coupling& coupling ) const 
{
  std::vector<std::pair<Particle, TotalCoupling>> rt;
  rt.emplace_back( coupling.particle(), TotalCoupling(coupling) );
  auto canDecay  = [&](const auto& particle){ return this->hasDecay(particle->name()); };
  auto canExpand = [&](const auto& coupling)
  {
    auto fs = coupling.getFinalStateParticles();
    return std::any_of( fs.begin(), fs.end(), canDecay);
  };

  do {
    std::vector<std::pair<Particle, TotalCoupling>> newTmpParticles;
    for(const auto& [particle, coupling] : rt )
    {
      if(!canExpand(particle)){ newTmpParticles.emplace_back(particle, coupling);  continue; }
      std::string nameToExpand = particle.uniqueString();
      auto fs                  = particle.getFinalStateParticles();
      auto it                  = *std::find_if( fs.begin(), fs.end(), canDecay );
      auto expandedRules       = rulesForDecay( it->name() );
      for ( auto& subTree : expandedRules )
      {
        auto expanded_amplitude = replaceAll( nameToExpand, it->name(), subTree.name() );
        newTmpParticles.emplace_back( Particle( expanded_amplitude), TotalCoupling( coupling, subTree) );
      }
    }
    rt = newTmpParticles;
  } while ( std::any_of( std::begin(rt), std::end(rt), [&](auto& it){ return canExpand(it.first);} ) );
  return rt;
}


void AmplitudeRules::add_rule( const Particle& p , double coupling ){ 
  m_rules[p.name()].emplace_back( p, coupling ); 
}

Coupling::Coupling(const Particle& particle, double f) :       
  m_name(particle.decayDescriptor()), 
  m_expr(new MinuitExpression(particle.decayDescriptor(), f)), 
  m_particle(particle){}
