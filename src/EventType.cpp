#include <math.h>
#include <stddef.h>
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

#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Units.h"
#include "AmpGen/Event.h"

using namespace AmpGen;
std::string convertTeXtoROOT(std::string input);

EventType::EventType( const std::vector<std::string>& particleNames, const bool& isTD ) : m_timeDependent( isTD )
{
  if ( particleNames.size() < 3 ) { // Mother plus two daughters minimum required
    ERROR( "Not enough particles in event type: " << particleNames[0] << " size =  " << particleNames.size() );
    throw std::runtime_error( "Not enough particles listed in particle names! Was it defined?" );
  }
  m_mother           = particleNames.at(0);
  for ( unsigned int i  = 1; i < particleNames.size(); ++i ) m_particleNames.push_back( particleNames[i] );
  auto motherProperties = ParticlePropertiesList::get( m_mother );
  if ( motherProperties != nullptr )
    m_motherMass = motherProperties->mass();
  else {
    ERROR( "Particle not found: " << m_mother );
    return;
  }
  m_alt_part_names = NamedParameter<bool>("EventType::AlternativeParicleNames", false );
  for ( auto& particle : m_particleNames ) {
    auto prop = ParticlePropertiesList::get( particle );
    if ( prop != nullptr )
      m_particleMasses.push_back( prop->mass() );
    else {
      ERROR( "Particle not found: " << particle );
      return;
    }
    if(m_alt_part_names)
      m_particleNamesPickled.push_back( replaceAll( replaceAll( particle, "+", "p" ), "-", "m" ) );
    else
      m_particleNamesPickled.push_back( replaceAll( replaceAll( particle, "+", "~" ), "-", "#" ) );
  }
  DEBUG( m_mother << " = " << m_motherMass << " -> " );
  for ( unsigned int i = 0; i < m_particleNames.size(); ++i ) {
    DEBUG( m_particleNames[i] << " = " << m_particleNamesPickled[i] << " = " << m_particleMasses[i] );
  }
  auto dimOfParticle = [](auto& name){
    return name == "gamma0" ? 2 : ParticlePropertiesList::get(name)->twoSpin() + 1;
  };
  m_dim.first = dimOfParticle(m_mother);
  m_dim.second = 1;
  for( auto& p : m_particleNames ) m_dim.second *= dimOfParticle(p);
}

std::map<std::string, size_t> EventType::getEventFormat( const bool& outputNames ) const
{
  std::map<std::string, size_t> returnValue;
  bool include_energy = NamedParameter<bool>("EventType::IncludeEnergy", true );
  size_t s = include_energy ? 4 : 3;
  for ( unsigned int ip = 0; ip < size(); ++ip ) {
    std::string parsed_name;
    if(m_alt_part_names)
      //check if there are multiple identical particles
      if(std::count(m_particleNamesPickled.begin(),m_particleNamesPickled.end(),m_particleNamesPickled[ip]) > 1)
        //if yes, append an index
        parsed_name = m_particleNamesPickled[ip] +
                      std::to_string(std::count(m_particleNamesPickled.begin(),m_particleNamesPickled.begin()+ip, m_particleNamesPickled[ip]));
      else // just take the already chosen name
        parsed_name = m_particleNamesPickled[ip];
    else
      parsed_name = "_" + std::to_string( ip + 1 ) + "_" + m_particleNamesPickled[ip];
    std::string stub = outputNames ? parsed_name : std::to_string( ip );
    if( include_energy ) returnValue[stub + "_E"]  = s * ip + 3;
    returnValue[stub + "_Px"] = s * ip + 0;
    returnValue[stub + "_Py"] = s * ip + 1;
    returnValue[stub + "_Pz"] = s * ip + 2;
  }
  if ( m_timeDependent ) returnValue[m_mother + "_decayTime"] = 4 * size();
  for( auto& extend : m_eventTypeExtensions ) returnValue[extend] = returnValue.size();
  return returnValue;
}
void EventType::extendEventType( const std::string& branch )
{
  m_eventTypeExtensions.push_back(branch);
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
std::pair<size_t, size_t> EventType::count(const size_t& index) const
{
  if( index >= size() ){
    ERROR("Looking for matching particles to index = " << index << " > size of eventType");
    return std::pair<size_t, size_t>(0, 0);
  }
  std::pair<size_t,size_t> rt(0,0);
  for( size_t j = 0 ; j < size(); ++j ){
    if( EventType::operator[](j) == EventType::operator[](index) ){
      rt.second++;
      if( j < index ) rt.first++;
    }
  }
  return rt;
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
    thing += label(x, isRoot) +" ";
  }
  return thing;
}

EventType EventType::conj( const bool& headOnly, const bool& dontConjHead ) const
{
  std::vector<std::string> type;
  type.push_back( dontConjHead ? m_mother : ParticlePropertiesList::get( m_mother )->anti().name() );
  std::transform( m_particleNames.begin(), m_particleNames.end(), std::back_inserter(type),
      [&](auto& x){ return headOnly ? x : ParticlePropertiesList::get(x)->anti().name() ; } );
  return EventType( type );
}

std::vector<Projection> EventType::defaultProjections(const size_t& nBins) const
{
  std::string defaultObservable  = NamedParameter<std::string>( "EventType::Observable", "mass2");
  std::vector<Projection> projections;
  for ( size_t r = 2; r < size(); ++r ) { /// loop over sizes ///
    std::vector<std::vector<size_t>>  combR = nCr( size(), r );
    std::transform( combR.begin(), combR.end(), std::back_inserter(projections),
      [&](auto& index){ return this->projection(nBins, index, defaultObservable ); } );
  }
  return projections;
}

Projection EventType::projection(const size_t& nBins, const std::vector<size_t>& indices, const std::string& observable) const
{
  bool useRootLabelling = NamedParameter<bool>("EventType::UseRootTEX", false );
  auto mm               = minmax(indices, true);
  std::string gevcccc   = useRootLabelling ? "GeV^{2}/c^{4}" : "\\mathrm{GeV}^{2}/c^{4}";
  std::string gevcc     = useRootLabelling ? "GeV/c^{2}"     : "\\mathrm{GeV}/c^{2}";
  if( observable == "mass2" )
    return Projection( [indices]( const Event& evt ) { return evt.s( indices ); },
        "s" + vectorToString( indices ),
        "s_{" + label( indices ) + "}", nBins,
        ( mm.first - 0.05 ) ,  ( mm.second + 0.05 ) , gevcccc );
  else if( observable == "mass" ){
    return Projection( [indices]( const Event& evt ) { return sqrt( evt.s( indices ) ); },
        "m" + vectorToString( indices ),
        "m_{" + label( indices ) + "}", nBins,
        mm.first > 0.05 ? sqrt(mm.first - 0.05) :0 ,  sqrt( mm.second + 0.05 ) , gevcc );
  }
  return Projection();
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
  os << type.mother() << " → [" << vectorToString(type.finalStates() , ", ") << "]";
  return os;
}

size_t EventType::dof() const { return 3 * size() - 7; }

std::function<void( Event& )> EventType::symmetriser() const
{
  std::map<std::string, std::vector<size_t>> particleOrdering;
  for ( size_t i = 0; i < m_particleNames.size(); ++i )
    particleOrdering[m_particleNames[i]].push_back( i );
  std::vector<std::vector<size_t>> shuffles;
  for ( auto& im : particleOrdering )
    if ( im.second.size() != 1 ) shuffles.push_back( im.second );

  int seed      = NamedParameter<unsigned int>( "EventType::SymmetriserSeed", 12 );
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
  return std::any_of( m_particleNames.begin(), m_particleNames.end(), [&name](auto& it) { return it == name ; } );
}

bool EventType::isTimeDependent() const { return m_timeDependent; }

size_t EventType::eventSize() const { return 4 * size() + m_timeDependent; }

std::pair<size_t, size_t> EventType::dim() const { return m_dim; }

std::string convertTeXtoROOT( std::string input )
{
  input = replaceAll( input, "\\mathrm{K}", "K" );
  input = replaceAll( input, "\\", "#" );
  input = replaceAll( input, "#xspace", "" );
  input = replaceAll( input, "#kern0.2em#overline{#kern-0.2em", "#bar{" );
  input = replaceAll( input, "^*", "^{*}" );
  return input;
}

std::string EventType::decayDescriptor() const
{
  return mother()+"{" + vectorToString(m_particleNames,",") +"}" ;
}
