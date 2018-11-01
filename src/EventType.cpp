#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

std::map<std::string, size_t> EventType::getEventFormat( const bool& outputNames ) const
{
  std::map<std::string, size_t> returnValue;

  for ( unsigned int ip = 0; ip < size(); ++ip ) {
    std::string stub =
        outputNames ? "_" + std::to_string( ip + 1 ) + "_" + m_particleNamesPickled[ip] : std::to_string( ip );
    returnValue[stub + "_E"]  = 4 * ip + 3;
    returnValue[stub + "_Px"] = 4 * ip + 0;
    returnValue[stub + "_Py"] = 4 * ip + 1;
    returnValue[stub + "_Pz"] = 4 * ip + 2;
  }
  if ( m_timeDependent ) returnValue[m_mother + "_ctau"] = 4 * size();
  return returnValue;
}

std::pair<double, double> EventType::minmax( const std::vector<size_t>& indices, bool isGeV ) const
{

  std::vector<size_t> ivec( size() );
  std::iota( ivec.begin(), ivec.end(), 0 );
  double min( 0 );
  double max( motherMass() );
  for ( auto& x : indices ) min += mass( x );
  for ( auto& x : ivec )
    if ( std::find( indices.begin(), indices.end(), x ) == indices.end() ) max -= mass( x );
  return std::pair<double, double>( min * min / GeV, max * max / GeV );
}

std::vector<std::vector<unsigned int>> EventType::getBosePairs() const
{
  std::map<std::string, std::vector<unsigned int>> particleOrdering;
  for ( unsigned int i = 0; i < m_particleNames.size(); ++i ) particleOrdering[m_particleNames[i]].push_back( i );
  std::vector<std::vector<unsigned int>> orderings;
  for ( auto& im : particleOrdering )
    if ( im.second.size() != 1 ) orderings.push_back( im.second );
  return orderings;
}
std::vector<std::string> EventType::finalStates() const { return m_particleNames; }
std::vector<double> EventType::masses() const { return m_particleMasses; }

size_t EventType::size() const { return m_particleNames.size(); }
double EventType::mass( const size_t& index ) const { return m_particleMasses[index]; }
double EventType::motherMass() const { return m_motherMass; }

std::string EventType::mother() const { return m_mother; }
std::string EventType::operator[]( const size_t& index ) const { return m_particleNames[index]; }
std::string EventType::label( const size_t& index, bool isRoot ) const
{
  const std::string label = ParticlePropertiesList::get( m_particleNames[index] )->label();
  return isRoot ? convertTeXtoROOT( label ) : label;
}

std::string EventType::label( const std::vector<size_t>& index, bool isRoot ) const
{
  std::string thing = "";
  for ( auto& x : index ) {
    thing += label( x, isRoot );
  }
  return thing;
}

EventType::EventType( const std::vector<std::string>& particleNames, const bool& isTD ) : m_timeDependent( isTD )
{
  
  if ( particleNames.size() < 3 ) { // Mother plus two daughters minimum required
    ERROR( "Not enough particles in event type: " << particleNames[0] << " size =  " << particleNames.size() );
    throw std::runtime_error( "Not enough particles listed in particle names! Was it defined?" );
  }
  m_mother           = particleNames.at( 0 );
  for ( unsigned int i  = 1; i < particleNames.size(); ++i ) m_particleNames.push_back( particleNames[i] );
  auto motherProperties = ParticlePropertiesList::get( m_mother );
  if ( motherProperties != nullptr )
    m_motherMass = motherProperties->mass();
  else {
    ERROR( "Particle not found: " << m_mother );
    return;
  }
  for ( auto& particle : m_particleNames ) {
    auto prop = ParticlePropertiesList::get( particle );
    if ( prop != nullptr )
      m_particleMasses.push_back( prop->mass() );
    else {
      ERROR( "Particle not found: " << particle );
      return;
    }
    m_particleNamesPickled.push_back( replaceAll( replaceAll( particle, "+", "~" ), "-", "#" ) );
  }
  DEBUG( m_mother << " = " << m_motherMass << " -> " );
  for ( unsigned int i = 0; i < m_particleNames.size(); ++i ) {
    DEBUG( m_particleNames[i] << " = " << m_particleNamesPickled[i] << " = " << m_particleMasses[i] );
  }
}

EventType EventType::conj( const bool& headOnly, const bool& dontConjHead ) const
{
  std::vector<std::string> type;
  type.push_back( dontConjHead ? m_mother : ParticlePropertiesList::get( m_mother )->anti().name() );

  for ( auto& x : m_particleNames ) type.push_back( headOnly ? x : ParticlePropertiesList::get( x )->anti().name() );
  return EventType( type );
}

std::vector<Projection> EventType::defaultProjections( const unsigned int& nBins ) const
{
  bool          useRootLabelling = NamedParameter<bool>( "EventType::UseRootTEX", false );
  std::string unitName           = NamedParameter<std::string>( "EventType::Units", "MeV" );
  std::string defaultObservable  = NamedParameter<std::string>( "EventType::Observable", "mass2");

  std::vector<Projection> axes;
  std::vector<std::vector<size_t>> permutations;
  std::string gevcccc            = useRootLabelling ? "GeV^{2}/c^{4}" : "\\mathrm{GeV}^{2}/c^{4}";
  std::string gevcc              = useRootLabelling ? "GeV/c^{2}" : "\\mathrm{GeV}/c^{2}";
  double units                   = 1;
  for ( size_t r = 2; r < size(); ++r ) { /// loop over sizes ///
    std::vector<std::vector<size_t>>  combR = nCr( size(), r );
    for ( auto& indices : combR ) {
      auto mm = minmax( indices, true );
      if( defaultObservable == "mass2" )
        axes.emplace_back( [indices, units]( const Event& evt ) { return evt.s( indices ); },
                         "s" + vectorToString( indices ),
                         "s_{" + label( indices, useRootLabelling ) + "}\\, \\left[" + gevcccc + "\\right]", nBins,
                         ( mm.first - 0.05 ) ,  ( mm.second + 0.05 ) , gevcccc );
      else if( defaultObservable == "mass" )
        axes.emplace_back( [indices, units]( const Event& evt ) { return sqrt( evt.s( indices ) / ( units * units ) ); },
                         "m" + vectorToString( indices ),
                         "m_{" + label( indices, useRootLabelling ) + "}\\, \\left[" + gevcc + "\\right]", nBins,
                         sqrt( mm.first - 0.05 ) ,  sqrt( mm.second + 0.05 ) , gevcc );
    }
  }
  return axes;
}

bool EventType::operator==( const EventType& other ) const
{
  if ( m_mother != other.mother() || size() != other.size() ) return false;

  for ( size_t i = 0; i < m_particleNames.size(); ++i ) {
    if ( m_particleNames[i] != other[i] ) return false;
  }
  return true;
}

std::ostream& AmpGen::operator<<( std::ostream& os, const EventType& type )
{
  os << type.mother() << " → [";
  auto fs = type.finalStates();
  for ( unsigned int i = 0; i < fs.size(); ++i ) os << fs[i] << ( i == fs.size() - 1 ? "]" : ", " );
  return os;
}

size_t EventType::dof() const { return 3 * size() - 7; }

std::function<void( Event& )> EventType::symmetriser() const
{
  auto shuffles = getBosePairs();
  int seed      = NamedParameter<unsigned int>( "EventType::symmetriser::seed", 12 );
  std::mt19937 rng( seed );
  for ( auto& shuffle : shuffles ) {
    std::string shuffle_string = "";
    for ( auto& s : shuffle ) shuffle_string += std::to_string( s ) + " ";
    DEBUG( "Shuffle = " << shuffle_string );
  }
  auto fcn = [shuffles, rng]( auto& event ) mutable -> void {
    for ( auto shuffled : shuffles ) {
      for ( unsigned int index = 0; index < shuffled.size(); ++index ) {
        unsigned int j = std::uniform_int_distribution<int>( 0, index )( rng );
        if ( index == j ) continue;
        std::swap( shuffled[index], shuffled[j] );
        event.swap( shuffled[index], shuffled[j] );
      }
    }
  };
  return fcn;
}

bool EventType::has( const std::string& name ) const
{
  for ( auto& it : m_particleNames )
    if ( it == name ) return true;
  return false;
}

bool EventType::isTimeDependent() const { return m_timeDependent; }

size_t EventType::eventSize() const { return 4 * size() + m_timeDependent; }

std::pair< size_t, size_t > EventType::dim() const 
{
  size_t it = ParticlePropertiesList::get(m_mother)->twoSpin() + 1;
  size_t ft = 1;
  for( auto& p : m_particleNames ) 
    ft *= ParticlePropertiesList::get( p )->twoSpin() + 1; 
  return {it,ft};
}
