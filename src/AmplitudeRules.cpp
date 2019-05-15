#include "AmpGen/AmplitudeRules.h"

#include <memory.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

AmplitudeRule::AmplitudeRule( const std::string& reName, 
                              const std::map<std::string, MinuitParameter*>& mapping )
{
  auto tokens = split( reName, '_' );
  if ( tokens.size() == 3 ) {
    m_prefix = tokens[0];
    m_name   = tokens[1];
    auto ire = mapping.find( reName );
    auto iim = mapping.find( tokens[0] + "_" + tokens[1] + "_Im" );
    if ( iim == mapping.end() ) {
      ERROR( "Well-formed coupling not identified for:" << reName );
      return;
    }
    m_re     = ire->second;
    m_im     = iim->second;
    m_isGood = true;
    
  } else if ( tokens.size() == 2 ) {
    m_prefix = "";
    m_name   = tokens[0];
    auto ire = mapping.find( reName );
    auto iim = mapping.find( tokens[0] + "_Im" );
    if ( iim == mapping.end() ) {
      ERROR( "Well-formed coupling not identified for:" << reName );
      return;
    }
    m_re     = ire->second;
    m_im     = iim->second;
    m_isGood = true;
  } else {
    ERROR( "Too many tokens! " );
    m_isGood = false;
  }
  if ( m_isGood ) {
    size_t pos = find_next_of( m_name, {"[", "{"} );
    if ( pos == std::string::npos ) 
    {
      ERROR( "Does not seem to be well formed decay descriptor [" << reName << "]" );
      m_isGood = false;
    }
  }
  m_particle = Particle(m_name);
}

AmplitudeRules::AmplitudeRules( const MinuitParameterSet& mps )
{
  for ( auto& o : mps.const_map() ) {
    if ( o.first.find( "_Re" ) == std::string::npos ) continue;
    AmplitudeRule pAmp( o.first, mps.const_map() );
    if ( pAmp.m_isGood ) m_rules[pAmp.m_particle.name()].push_back( pAmp );
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
  std::vector<std::string> particleNames;
  particleNames.push_back( particle.name() );
  std::vector<std::shared_ptr<Particle>> fs = particle.getFinalStateParticles();
  std::stable_sort( fs.begin(), fs.end(), []( auto& A, auto& B ) { return *A < *B; } );
  for( auto& f : fs ) particleNames.push_back( f->name() );
  return EventType( particleNames );
}

CouplingConstant::CouplingConstant(const AmplitudeRule& pA)
{
  couplings.emplace_back( std::make_pair(pA.m_re,pA.m_im) );  
  std::string cartOrPolar = NamedParameter<std::string>("CouplingConstant::Coordinates" ,"cartesian");
  std::string degOrRad    = NamedParameter<std::string>("CouplingConstant::AngularUnits","rad");
  if( cartOrPolar == "polar" ){
    isCartesian = false; 
  }
  else if ( cartOrPolar != "cartesian" ){
    FATAL("Coordinates for coupling constants must be either cartesian or polar");
  } 
  if ( degOrRad == "deg") sf = M_PI / 180; 
  else if ( degOrRad != "rad"){
    FATAL("CouplingConstant::AngularUnits must be either rad or deg");
  } 
}

std::complex<double> CouplingConstant::operator()() const
{
  std::complex<double> F( 1, 0 );
  if ( isCartesian )
    for( auto& p : couplings ) F *= complex_t( p.first->mean() , p.second->mean() ); 
  else
    for( auto& p : couplings ) F *= p.first->mean() * complex_t( cos( sf * p.second->mean() ), sin( sf * p.second->mean() ) );
  return F;
}

Expression CouplingConstant::to_expression() const
{
  Expression J = Constant(0,1);
  if ( isCartesian ) {
    Expression total = 1;
    for ( auto& p : couplings ) {
      total = total * ( Parameter( p.first->name() ) + J * Parameter( p.second->name() ) );
    }
    return total;
  } else {
    Expression angle = 0;
    Expression amp   = 1;
    for ( auto& p : couplings ) {
      angle = angle + Parameter( p.second->name() );
      amp   = amp * Parameter( p.first->name() );
    }
    return amp * ( Cos( sf * angle ) + J * Sin( sf * angle ) );
  }
}

void CouplingConstant::print() const
{
  INFO( couplings[0].first->name() << " (" << isCartesian << ") = " << ( *this )() );
  if ( isCartesian )
    for ( auto& coupling : couplings ) INFO( coupling.first->name() + " i " + coupling.second->name() );
  else
    for ( auto& coupling : couplings ) INFO( coupling.first->name() << " x exp(i" <<  coupling.second->name() );
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
  for( auto& c : couplings ) 
    if( c.first->iFixInit() == 0 or c.second->iFixInit() == 0 ) return false; 
  return true;
}

bool CouplingConstant::contains( const std::string& label ) const 
{
  for ( auto& reAndIm : couplings )
    if ( reAndIm.first->name().find( label ) != std::string::npos ) return true; 
  return false; 
}

void CouplingConstant::changeSign() 
{
  auto& top_coupling = *couplings.begin();
  if( isCartesian ) top_coupling.first->setCurrentFitVal( top_coupling.first->mean() * -1. );
}

