#include <cmath>
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
#include "AmpGen/OptionsParser.h"
#include "AmpGen/PhaseSpace.h"
using namespace AmpGen;

EventType::EventType( const std::vector<std::string>& particleNames, const bool& isTD ) : m_timeDependent( isTD )
{
  if ( OptionsParser::printHelp() ) return; 

  if ( particleNames.size() < 3 ) { // Mother plus two daughters minimum required
    ERROR( "Not enough particles in event type: " << particleNames[0] << " size =  " << particleNames.size() );
    throw std::runtime_error( "Not enough particles listed in particle names! Was it defined?" );
  }
  m_mother           = particleNames.at(0);
  m_alt_part_names = NamedParameter<bool>("EventType::AlternativeParticleNames", false );
  
  auto motherProperties = ParticlePropertiesList::get( m_mother );
  if ( motherProperties != nullptr )
    m_motherMass = motherProperties->mass();
  else {
    ERROR( "Particle not found: " << m_mother );
    return;
  }
  
  m_ignore.resize(particleNames.size()-1,false);
  for ( auto it = particleNames.begin()+1; it != particleNames.end(); ++it )
  { 
    if( *it->begin() == '{' && *it->rbegin() == '}' )
    {
      m_particleNames.push_back( it->substr(1, it->size() -2 ) );
    }
    else m_particleNames.push_back( *it );

    auto prop = ParticlePropertiesList::get( *m_particleNames.rbegin() );
    if ( prop != nullptr )
      m_particleMasses.push_back( prop->mass() );
    else {
      FATAL( "Particle not found: " << *m_particleNames.rbegin() );
      return;
    }
    if(m_alt_part_names)
      m_particleNamesPickled.push_back( replaceAll( *m_particleNames.rbegin(), {{"+","p"},{"-","m"},{"(",""},{")",""}}));
    else
      m_particleNamesPickled.push_back( replaceAll( *m_particleNames.rbegin(), {{"+","~"},{"-","#"},{"(",""},{")",""}}));
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

std::map<std::string, unsigned> EventType::getEventFormat( const bool& outputNames ) const
{
  std::map<std::string, unsigned> returnValue;
  bool include_energy = NamedParameter<bool>("EventType::IncludeEnergy", true );
  unsigned s = include_energy ? 4 : 3;
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
std::pair<double, double> EventType::minmax( const std::vector<unsigned>& indices) const
{
  std::vector<unsigned> ivec( size() );
  std::iota( ivec.begin(), ivec.end(), 0 );
  double min = 0 ; 
  for( auto& i : indices ) min += mass(i);
  double max = motherMass();
  for ( auto& x : ivec )
    if ( std::find( indices.begin(), indices.end(), x ) == indices.end() ) max -= mass( x );
  return std::pair<double, double>(min, max);
}
std::pair<unsigned, unsigned> EventType::count(const unsigned& index) const
{
  if( index >= size() ){
    ERROR("Looking for matching particles to index = " << index << " > size of eventType");
    return std::pair<unsigned, unsigned>(0, 0);
  }
  std::pair<unsigned,unsigned> rt(0,0);
  for( unsigned j = 0 ; j < size(); ++j ){
    if( EventType::operator[](j) == EventType::operator[](index) ){
      rt.second++;
      if( j < index ) rt.first++;
    }
  }
  return rt;
}

std::vector<std::string> EventType::finalStates() const { return m_particleNames; }
std::vector<double> EventType::masses() const { return m_particleMasses; }

unsigned EventType::size() const { return m_particleNames.size(); }
double EventType::mass( const unsigned& index ) const { return m_particleMasses[index]; }
double EventType::motherMass() const { return m_motherMass; }

std::string EventType::mother() const { return m_mother; }
std::string EventType::operator[]( const unsigned& index ) const { return m_particleNames[index]; }
std::string EventType::label( const unsigned& index, bool isRoot ) const
{
  return ParticlePropertiesList::get( m_particleNames[index] )->label();
}

std::string EventType::label( const std::vector<unsigned>& index, bool isRoot ) const
{
  std::string thing = "";
  for ( auto& x : index ) thing += label(x, isRoot) +" ";
  return thing;
}

EventType EventType::conj( const bool& headOnly, const bool& dontConjHead ) const
{
  std::vector<std::string> type;
  type.push_back( dontConjHead ? m_mother : ParticlePropertiesList::get( m_mother )->anti().name() );
  std::transform( m_particleNames.begin(), m_particleNames.end(), std::back_inserter(type),
      [&](auto& x){ return headOnly ? x : ParticlePropertiesList::get(x)->anti().name() ; } );
  return EventType( type, m_timeDependent );
}

std::vector<Projection> EventType::defaultProjections(const unsigned& nBins) const
{
  std::string defaultObservable  = NamedParameter<std::string>( "EventType::Observable", "mass2");
  std::vector<Projection> projections;
  for ( unsigned r = 2; r < size(); ++r ) { /// loop over sizes ///
    std::vector<std::vector<unsigned>>  combR = nCr( size(), r );
    std::transform( combR.begin(), combR.end(), std::back_inserter(projections),
      [&](auto& index){ return this->projection(nBins, index, defaultObservable ); } );
  }
  return projections;
}

Projection EventType::projection(const unsigned& nBins, const std::vector<unsigned>& indices, const std::string& observable) const
{
  bool useRootLabelling = NamedParameter<bool>("EventType::UseRootTEX", false );
  auto mm               = minmax(indices);
  std::string gevcccc   = useRootLabelling ? "GeV^{2}/c^{4}" : "\\mathrm{GeV}^{2}/c^{4}";
  std::string gevcc     = useRootLabelling ? "GeV/c^{2}"     : "\\mathrm{GeV}/c^{2}";
  if( observable == "mass2" )
    return Projection( [indices]( const Event& evt ) { return evt.s( indices ); },
        "s" + vectorToString( indices ),
        "s_{" + label( indices ) + "}", nBins,
        ( mm.first * mm.first - 0.05 ) ,  ( mm.second * mm.second + 0.05 ) , gevcccc );
  else if( observable == "mass" ){
    return Projection( [indices]( const Event& evt ) { return sqrt( evt.s( indices ) ); },
        "m" + vectorToString( indices ),
        "m_{" + label( indices ) + "}", nBins,
        mm.first > 0.05 ? mm.first - 0.05 :0 ,  mm.second + 0.05, gevcc );
  }
  return Projection();
}

bool EventType::operator==( const EventType& other ) const
{
  if ( m_mother != other.mother() || size() != other.size() ) return false;
  for ( unsigned i = 0; i < m_particleNames.size(); ++i ) {
    if ( m_particleNames[i] != other[i] ) return false;
  }
  return true;
}

std::ostream& AmpGen::operator<<( std::ostream& os, const EventType& type )
{
  os << type.mother() << " → [" << vectorToString(type.finalStates() , ", ") << "]";
  return os;
}

unsigned EventType::dof() const { return 3 * size() - 7; }

std::function<void( Event& )> EventType::symmetriser() const
{
  std::map<std::string, std::vector<unsigned>> particleOrdering;
  for ( unsigned i = 0; i < m_particleNames.size(); ++i )
    particleOrdering[m_particleNames[i]].push_back( i );
  std::vector<std::vector<unsigned>> shuffles;
  for ( auto& im : particleOrdering )
    if ( im.second.size() != 1 ) shuffles.push_back( im.second );

  int seed      = NamedParameter<unsigned int>( "EventType::SymmetriserSeed", 12 );
  std::mt19937 rng( seed );
  for ( auto& shuffle : shuffles ) {
    std::string shuffle_string = "";
    for ( auto& s : shuffle ) shuffle_string += std::to_string( s ) + " ";
    DEBUG( "Shuffle = " << shuffle_string );
  }
  return [shuffles, rng]( auto& event ) mutable -> void {
    for ( auto shuffled : shuffles ) {
      for ( unsigned int index = 0; index < shuffled.size(); ++index ) {
        unsigned int j = std::uniform_int_distribution<int>( 0, index )( rng );
        if ( index == j ) continue;
        std::swap( shuffled[index], shuffled[j] );
        event.swap( shuffled[index], shuffled[j] );
      }
    }
  };
}

std::function<bool(Event&, const std::vector<int>&)> EventType::automaticOrdering() const 
{
  std::vector<int> ids;
  for( unsigned i = 0 ; i != m_particleNames.size(); ++i ) ids.push_back( ParticleProperties::get(m_particleNames[i])->pdgID() );
  auto matches = [](const auto& c1, const auto& c2 , unsigned sgn = +1)
  {
    std::vector<bool> used( c1.size(), false );
    for(unsigned i = 0; i != c1.size(); ++i )
    {
      for( unsigned j = 0; j != c2.size(); ++j )
      {
        if( c1[i] == sgn * c2[j] && ! used[j] ) used[j] = true;
      }
    }
    return std::all_of( std::begin(used), std::end(used), [](auto b) { return b; }  ) ;
  };

  return [ids, matches](auto& event, const auto& actual_ids) -> bool {
    std::vector<unsigned> new_addresses( ids.size(), 999 ); 
    int sgn = +1; 
    if( matches(ids, actual_ids ) ) sgn = +1;
    else if( matches(ids, actual_ids, -1 ) ) sgn = -1;
    else { ERROR("Ids: " << vectorToString(actual_ids, " ") << " do not match either particle or antiparticle ["<< vectorToString(ids, " ") << "]" );
      return false; 
    }
    
    for( unsigned i = 0 ; i != ids.size(); ++i )
    {
      for( unsigned j = 0 ; j != actual_ids.size(); ++j )
      {
        if( actual_ids[j] ==  sgn * ids[i] && new_addresses[j]==999 ){ new_addresses[j] = i; break ; }
      }
    }
    event.reorder( new_addresses );
    return true; 
  };
}


bool EventType::has( const std::string& name ) const
{
  return std::any_of( m_particleNames.begin(), m_particleNames.end(), [&name](auto& it) { return it == name ; } );
}

bool EventType::isTimeDependent() const { return m_timeDependent; }

unsigned EventType::eventSize() const { return 4 * size() + m_timeDependent; }

std::pair<unsigned, unsigned> EventType::dim() const { return m_dim; }

std::string EventType::decayDescriptor() const
{
  return mother()+"{" + vectorToString(m_particleNames,",") +"}" ;
}

Event EventType::makeEvent() const 
{
  return PhaseSpace( *this ).makeEvent(); 
}

std::string EventType::constructor_string() const 
{
  return mother() + " " + vectorToString(m_particleNames, " "); 
}

extern "C" unsigned AmpGen::python__EventType__dim(const char* eventType){
  
  EventType type( split( std::string(eventType),' ') );
  auto dim = type.dim();
  return (dim.first << 16) + dim.second;   
}
