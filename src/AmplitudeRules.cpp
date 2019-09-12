#include "AmpGen/AmplitudeRules.h"

#include <memory.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <numeric>

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;
using namespace std::complex_literals; 

AmplitudeRule::AmplitudeRule(MinuitParameter* re, MinuitParameter* im ) :
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
}

AmplitudeRules::AmplitudeRules( const MinuitParameterSet& mps )
{
  for ( auto& it_re : mps ) {
    if ( it_re->name().find("_Re") == std::string::npos ) continue;
    auto it_im = mps.find(replaceAll( it_re->name(), "_Re","_Im") );
    if( it_im == nullptr ){
      ERROR("Cannot find matching imaginary part / phase for: " <<  it_re->name() );
      continue; 
    }
    auto p = AmplitudeRule(it_re, it_im );
    m_rules[p.head()].emplace_back(p);
  }
}

CouplingConstant::CouplingConstant( const CouplingConstant& other, 
                    const AmplitudeRule& pA) : 
  couplings(other.couplings ),
  isCartesian(other.isCartesian),
  sf(other.sf)
{
  couplings.emplace_back( pA.m_re, pA.m_im );
}


bool AmplitudeRules::hasDecay(const std::string& head) 
{ 
  return m_rules.find(head) != m_rules.end(); 
}

std::vector<AmplitudeRule> AmplitudeRules::rulesForDecay(const std::string& head, const std::string&  prefix)
{
  if(!hasDecay(head)) return std::vector<AmplitudeRule>();
  if( prefix == "" )return m_rules[head];
  std::vector<AmplitudeRule> rt = m_rules[head];
  rt.erase( std::remove_if( std::begin(rt), std::end(rt), [&prefix](auto& p){ return p.prefix() != prefix; } ) );
  return rt;
}

std::map<std::string, std::vector<AmplitudeRule>> AmplitudeRules::rules() 
{ 
  return m_rules;
}

EventType AmplitudeRule::eventType() const
{
  Particle particle( m_name );
  std::vector<std::string> particleNames = { particle.name() };
  std::vector<std::shared_ptr<Particle>> fs = particle.getFinalStateParticles();
  std::stable_sort( fs.begin(), fs.end(), []( auto& A, auto& B ) { return *A < *B; } );
  std::transform( fs.begin(), fs.end(), std::back_inserter(particleNames), [](auto& p ) -> std::string { return p->name() ; } );
  return EventType( particleNames );
}

namespace AmpGen 
{
  make_enum(coordinateType, cartesian, polar)
  make_enum(angType, deg, rad)
}

CouplingConstant::CouplingConstant(const AmplitudeRule& pA)
{
  couplings.emplace_back(pA.m_re,pA.m_im);  
  coordinateType coord = NamedParameter<coordinateType>("CouplingConstant::Coordinates", coordinateType::cartesian);
  angType degOrRad     = NamedParameter<angType>("CouplingConstant::AngularUnits", angType::rad);
  if( coord == coordinateType::polar ) isCartesian = false; 
  else if ( coord != coordinateType::cartesian){
    FATAL("Coordinates for coupling constants must be either cartesian or polar");
  } 
  if ( degOrRad == angType::deg) sf = M_PI / 180; 
  else if ( degOrRad != angType::rad){
    FATAL("CouplingConstant::AngularUnits must be either rad or deg");
  } 
}

std::complex<double> CouplingConstant::operator()() const
{
  return isCartesian ?  
    std::accumulate( couplings.begin(), couplings.end(), complex_t(1,0), 
        [](auto& prod, auto& coupling){ return prod * complex_t( coupling.first->mean(), coupling.second->mean() ) ; } )
   : std::accumulate( couplings.begin(), couplings.end(), complex_t(1,0), 
        [this](auto& prod, auto& coupling){ return prod * coupling.first->mean() * exp( 1i* this->sf * coupling.second->mean() ) ; } );
}

Expression CouplingConstant::to_expression() const
{
  if ( isCartesian ) {
    return std::accumulate( couplings.begin(), couplings.end(), Expression(1), 
        [](auto& prod, auto& p){ return prod * ( Parameter(p.first->name()) + 1i*Parameter(p.second->name() ) ) ; } );
  } else {
    auto it = std::accumulate( couplings.begin(), couplings.end(), std::pair<Expression,Expression>(1,0),
        [](auto& prod, auto& p) -> std::pair<Expression, Expression> { 
        return std::make_pair( prod.first * Parameter( p.first->name() ), prod.second + Parameter( p.second->name() ) ); } );
    return it.first * fcn::exp( 1i*sf*it.second );
  }
}

void CouplingConstant::print() const
{
  INFO( couplings[0].first->name() << " (" << isCartesian << ") = " << ( *this )() );
  if ( isCartesian )
    for ( auto& coupling : couplings ) INFO( coupling.first->name() + " i " + coupling.second->name() );
  else
    for ( auto& coupling : couplings ) 
      INFO( coupling.first->name() << " x exp(i" <<  coupling.second->name() << ") = " << 
          coupling.first->mean() * exp( 1i * coupling.second->mean() * M_PI / 180. )   );
}

std::vector<std::pair<Particle, CouplingConstant>> AmplitudeRules::getMatchingRules(const EventType& type, const std::string& prefix )
{
  auto rules        = rulesForDecay( type.mother() );
  std::vector<std::pair<Particle, CouplingConstant>> rt; 
  for ( auto& rule : rules ) {
    if ( rule.prefix() != prefix ) continue;
    std::vector<std::pair<Particle, CouplingConstant>> tmpParticles;
    auto fs = type.finalStates();
    tmpParticles.emplace_back( Particle( rule.name(), fs ), CouplingConstant(rule) );
    do {
      std::vector<std::pair<Particle, CouplingConstant>> newTmpParticles;
      for ( auto& particleWithCouplingConstant : tmpParticles ) {
        auto protoParticle    = particleWithCouplingConstant.first;
        auto coupling         = particleWithCouplingConstant.second;
        auto protoFinalStates = protoParticle.getFinalStateParticles();
        if ( protoFinalStates.size() == type.size() ) {
          rt.emplace_back( particleWithCouplingConstant );
          continue;
        }
        std::string nameToExpand = protoParticle.uniqueString();
        for ( auto& ifs : protoFinalStates ) {
          auto expandedRules = rulesForDecay( ifs->name() ); /// get rules for decaying particle
          if ( expandedRules.size() == 0 ) continue;
          for ( auto& subTree : expandedRules ) {
            auto expanded_amplitude = replaceAll( nameToExpand, ifs->name(), subTree.name() );
            auto fs2                = type.finalStates();
            newTmpParticles.emplace_back( Particle( expanded_amplitude, fs2 ), CouplingConstant( coupling, subTree) );
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

bool CouplingConstant::isFixed() const
{
  return std::all_of( couplings.begin(), couplings.end(), 
      [](auto& c){ return c.first->isFixed() && c.second->isFixed() ; } );
}

bool CouplingConstant::contains( const std::string& label ) const 
{
  return std::any_of( couplings.begin(), couplings.end(), 
    [&label](auto& c){ return c.first->name().find(label) != std::string::npos ; } );
}

