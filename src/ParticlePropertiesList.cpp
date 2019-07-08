// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:04 GMT
#include "AmpGen/ParticlePropertiesList.h"

#include <stdlib.h>
#include <stdexcept>
#include <string>

#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

ParticlePropertiesList* ParticlePropertiesList::ptr = nullptr;

const ParticlePropertiesList* ParticlePropertiesList::getMe()
{
  if(!ptr) ptr = new ParticlePropertiesList();
  if(nullptr == ptr) FATAL("Couldn't get ParticlePropertiesList (i.e. myself)" );
  return ptr;
}

ParticlePropertiesList* ParticlePropertiesList::getMutable()
{
  if(!ptr) ptr = new ParticlePropertiesList();
  if(ptr == nullptr) FATAL( "Couldn't get ParticlePropertiesList (i.e. myself)" );
  return ptr;
}

const ParticleProperties* ParticlePropertiesList::get( const std::string& name, const bool& quiet )
{
  const ParticleProperties* props = getMe()->find( name, quiet );
  if ( nullptr == props ) return nullptr;
  return props;
}

const ParticleProperties* ParticlePropertiesList::get( const int& pid, const bool& quiet )
{
  const ParticleProperties* props = getMe()->find( pid, quiet );
  if ( nullptr == props ) return nullptr;
  return props;
}

const std::vector<std::string> ParticlePropertiesList::dirList() const
{
  std::vector<std::string> dirList;

  std::string AmpGenRoot( "." );
  char* AmpEnv = getenv( "AMPGENROOT" );
  if ( nullptr != AmpEnv ) {
    AmpGenRoot = AmpEnv;
  } else {
#ifdef AMPGENROOT_CMAKE
    AmpGenRoot = AMPGENROOT_CMAKE;
    INFO( "Using built-in AMPGENROOT (set env variable if incorrect): " << AmpGenRoot );
#else
    WARNING( "AMPGENROOT not set, may not be able to find pdg database" );
#endif
  }
  dirList.push_back( "" );
  dirList.push_back( AmpGenRoot + "/" );
  dirList.push_back( AmpGenRoot + "/options/" );
  dirList.push_back( "../" );
  return dirList;
}

ParticlePropertiesList::ParticlePropertiesList( const std::string& fname_in )
{
  auto dl = dirList();
  bool status = true; 
  status &= std::any_of( dl.begin(), dl.end(), [this](auto& d){ return this->readLatexLabels(d +"pdgID_to_latex.dat") ; } );
  status &= std::any_of( dl.begin(), dl.end(), [this](auto& d){ return this->readFile(d +"mass_width.csv") ; } );
  status &= std::any_of( dl.begin(), dl.end(), [this](auto& d){ return this->readFile(d +"MintDalitzSpecialParticles.csv") ; } );
  if( !status ){
    WARNING("Failed to load full PDG configuration, beware of unexpected behaviour");
  }
  makeMappings();
  m_quasiStableThreshold = NamedParameter<double>( "ParticleProperties::qsThreshold", KeV ); /// limit is 1 keV
}

bool ParticlePropertiesList::readLatexLabels( const std::string& name )
{
  if ( !fileExists( name ) ) return false;
  m_latexLabels.clear();
  processFile( name, [this]( auto& line ) {
    auto tokens                            = split( line, ' ' );
    this->m_latexLabels[stoi( tokens[0] )] = std::make_pair( tokens[1], tokens[2] );
  } );
  return true;
}

bool ParticlePropertiesList::readFile( const std::string& name )
{

  if ( !fileExists( name ) ) {
    DEBUG( "File not found: " << name );
    return false;
  }
  INFO( "Reading file: " << name );
  processFile( name, [this]( auto& line ) {
    if ( line[0] == '*' ) return;
    ParticleProperties P( line );
    if ( !P.isValid() ) {
      WARNING( line << " is not valid" );
      return;
    }
    auto label = m_latexLabels.find( P.pdgID() );
    if ( label != m_latexLabels.end() ) P.setLabel( label->second.first );
    m_theList.push_back( P );
    if ( P.hasDistinctAnti() ) {
      P.antiThis();
      if ( label != m_latexLabels.end() ) P.setLabel( label->second.second );
      m_theList.push_back( P );
    }
  } );
  return true;
}

void ParticlePropertiesList::makeMappings()
{
  for ( auto& it : m_theList ) {
    m_byName[it.name()] = &it;
    if ( it.pdgID() != 0 ) {
      auto found = m_byID.find( it.pdgID() );
      if ( found != m_byID.end() ) {
        WARNING( "pdgID " << it.pdgID() << " used twice, here: " << ( found->second )->name()
                           << ", and here: " << it.name() );
      }
    }
    m_byID[it.pdgID()] = &it;
  }
}
void ParticlePropertiesList::print( std::ostream& out ) const
{
  for ( auto& x : m_theList ) {
    x.print( out );
    out << "\n";
  }
}

const ParticleProperties* ParticlePropertiesList::find( const std::string& name, bool quiet ) const
{
  auto it = m_byName.find( name );
  if ( it != m_byName.end() ) return it->second;
  if ( !quiet ) {
    auto particleNames = ParticlePropertiesList::getMe()->getParticleNames();

    unsigned int minDistance = 9999;
    std::string suggestion   = "";
    for ( auto& particle : particleNames ) {
      unsigned int distance = editDistance( particle, name );
      if ( distance < minDistance ) {
        suggestion  = particle;
        minDistance = distance;
      }
    }
    ERROR( "Particle: " << name << " not in PDG. Did you mean " << suggestion << "?" );
  }
  return nullptr;
}
const ParticleProperties* ParticlePropertiesList::find( int id, bool quiet ) const
{
  auto it = m_byID.find( id );
  if ( it == m_byID.end() ) {
    if ( !quiet ) ERROR( "Particle with id: " << id << " not recognised" );
    return nullptr;
  }
  return it->second;
}

std::ostream& operator<<( std::ostream& out, const ParticlePropertiesList& ppl )
{
  ppl.print( out );
  return out;
}

std::vector<std::string> ParticlePropertiesList::getParticleNames() const
{
  std::vector<std::string> particleNames;
  std::transform( m_byName.begin(), m_byName.end(), std::back_inserter(particleNames), [](auto& p ) -> std::string { return p.first ; } );
  return particleNames;
}
std::vector<int> ParticlePropertiesList::getParticleIds() const
{
  std::vector<int> particleIds;
  std::transform( m_byID.begin(), m_byID.end(), std::back_inserter(particleIds), [](auto& p ) -> int { return p.first ; } );
  return particleIds;
}

void ParticlePropertiesList::makeAlias( const std::string& name, const std::string& alias )
{
  DEBUG( "Making alias: " << alias << " for " << name );
  auto it = find( name );
  if ( it == nullptr ) {
    ERROR( "Cannot find particle " << name << " with which to make alias" );
    return;
  }
  ParticleProperties* pp = new ParticleProperties( *it );
  pp->setName( alias );
  m_byName[alias] = pp;
}

double ParticlePropertiesList::quasiStableThreshold() const { return m_quasiStableThreshold; }
