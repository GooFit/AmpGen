#include "AmpGen/AmplitudeRules.h"

#include <memory.h>
#include <ext/alloc_traits.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

AmplitudeRules::AmplitudeRules( MinuitParameterSet& mps )
{

  for ( auto& o : mps.map() ) {
    if ( o.first.find( "_Re" ) == std::string::npos ) continue;
    AmplitudeRule pAmp( o.first, mps.map() );
    if ( pAmp.m_isGood ) {
      m_rules[pAmp.m_head].push_back( pAmp );
    }
  }
}

bool AmplitudeRules::hasDecay( const std::string& head ) { return m_rules.find( head ) != m_rules.end(); }

std::vector<AmplitudeRule> AmplitudeRules::rulesForDecay( const std::string& head )
{
  if ( !hasDecay( head ) ) return std::vector<AmplitudeRule>();
  return m_rules[head];
}

std::map<std::string, std::vector<AmplitudeRule>> AmplitudeRules::rules() { return m_rules; }

AmplitudeRule::AmplitudeRule( const std::string& reName, std::map<std::string, MinuitParameter*>& mapping )
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
    m_re     = mapping[reName];
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
    if ( pos != std::string::npos )
      m_head = m_name.substr( 0, pos );
    else {
      ERROR( "Does not seem to be well formed decay descriptor [" << reName << "]" );
      m_isGood = false;
    }
  }
}
EventType AmplitudeRule::eventType() const
{
  Particle particle( m_name );
  std::vector<std::string> particleNames;
  particleNames.push_back( particle.name() );
  std::vector<std::shared_ptr<Particle>> fs = particle.getFinalStateParticles();
  std::stable_sort( fs.begin(), fs.end(), []( auto& A, auto& B ) { return *A < *B; } );
  for ( auto& f : fs ) particleNames.push_back( f->name() );
  return EventType( particleNames );
}

Coupling AmplitudeRule::makeCoupling( bool isCartan )
{
  Coupling c;
  c.couplings.emplace_back( m_re, m_im );
  c.isCartesian = isCartan;
  return c;
}

std::complex<double> Coupling::operator()() const
{

  if ( isCartesian ) {
    std::complex<double> F( 1, 0 );
    for ( auto& p : couplings ) {
      // INFO( p.first << " " << p.second );
      F *= std::complex<double>( p.first->mean(), p.second->mean() );
    }
    return F;
  } else {
    std::complex<double> F( 1, 0 );
    for ( auto& p : couplings ) {
      F *= p.first->mean() * std::complex<double>( cos( p.second->mean() ), sin( p.second->mean() ) );
    }
    // INFO("Returning F = " << F );
    return F; // amp * std::complex<double>( cos(angle), sin(angle) );
  }
}

Expression Coupling::to_expression() const
{
  Expression J = Constant(0,1);
  if ( isCartesian ) {
    Expression total =1 ;
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
    return amp * ( Cos( angle ) + J * Sin( angle ) );
  }
}

void Coupling::print() const
{
  INFO( couplings[0].first->name() << " (" << isCartesian << ") = " << ( *this )() );
  if ( isCartesian )
    for ( auto& coupling : couplings ) INFO( coupling.first->name() + " i " + coupling.second->name() );
  else
    for ( auto& coupling : couplings ) INFO( coupling.first->name() << " x exp(i" << coupling.second->name() );
}

std::vector< std::pair<Particle, Coupling > > AmplitudeRules::getMatchingRules( 
  const AmpGen::EventType& type, const std::string& prefix, const bool& useCartesian ){

  auto rules        = rulesForDecay( type.mother() );
  std::vector< std::pair< Particle, Coupling > > rt; 

  for ( auto& p : rules ) {
    if ( p.prefix() != prefix ) continue;
    std::vector<std::pair<Particle, Coupling>> tmpParticles;
    auto fs = type.finalStates();
    tmpParticles.emplace_back( Particle( p.name(), fs ), p.makeCoupling( useCartesian ) );
    do {
      std::vector<std::pair<Particle, Coupling>> newTmpParticles;
      for ( auto& particleWithCoupling : tmpParticles ) {
        auto protoParticle    = particleWithCoupling.first;
        auto coupling         = particleWithCoupling.second;
        auto protoFinalStates = protoParticle.getFinalStateParticles();
        if ( protoFinalStates.size() == type.size() ) {
          rt.emplace_back( particleWithCoupling );
          //        addMatrixElement( particleWithCoupling );
          continue; /// this particle is fully expanded
        }
        std::string nameToExpand = protoParticle.uniqueString();
        for ( auto& ifs : protoFinalStates ) {
          auto expandedRules = rulesForDecay( ifs->name() ); /// get rules for decaying particle
          if ( expandedRules.size() == 0 ) continue;
          for ( auto& subTree : expandedRules ) {
            auto expanded_amplitude = replaceAll( nameToExpand, ifs->name(), subTree.name() );
            auto fs2                = type.finalStates();
            newTmpParticles.emplace_back( Particle( expanded_amplitude, fs2 ),
                                          Coupling( coupling, subTree, useCartesian ) );
          }
          break; // we should only break if there are rules to be expanded ...
        }
      }
      tmpParticles = newTmpParticles;
    } while ( tmpParticles.size() != 0 );
  }
  return rt;
}
