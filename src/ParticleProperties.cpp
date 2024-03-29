// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:04 GMT
#include "AmpGen/ParticleProperties.h"

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <cmath>

#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

void ParticleProperties::print( std::ostream& out ) const
{
  out << "Mass  " << mass() << " +" << mErrPlus() << " -" << mErrMinus() << "\nWidth " << width() << " +" << wErrPlus()
    << " -" << wErrMinus() << "\n I=" << I() << ", G=" << G() << "\n J=" << J() << ", C=" << C() << ", P=" << P()
    << "\n Q = " << charge() << "\n pdgID " << pdgID() << "\n name " << name()
    << "\n quark-content " << quarkContent() << "\n is its own antiparticle? "
    << ( hasDistinctAnti() ? "no" : "yes" ) << "\n radius " << radius() * GeV << " /GeV"
    << "\n";
}

std::vector<int> ParticleProperties::polarisations() const 
{
  if( twoSpin() == 0 ) return {0};          // scalar 
  if( isPhoton() )     return {1,-1};       // photon 
  if( twoSpin() == 1 ) return {1,-1};       // fermion
  if( twoSpin() == 2 ) return {1,0,-1};     // vector
  if( twoSpin() == 4 ) return {2,1,0,-1,-2};// tensor
  else { 
    WARNING("Particle with spin: " << twoSpin() << "/2" << " not implemented in initial/final state");
    return {0};
  }
}

int ParticleProperties::chargeFromString( const std::string& ch, bool& status ) const
{
  if ( ch == "+" ) return 1;
  if ( ch == "-" ) return -1;
  if ( ch == " " ) return 0;
  if ( ch == "0" )  return  0;
  if ( ch == "++" ) return  2;
  if ( ch == "--" ) return -2;
  status = 0;
  return 0;
}

std::string ParticleProperties::J() const 
{
  if( m_twoSpin == -1 ) return "?";
  if( m_twoSpin % 2 == 0 ) return std::to_string( int(m_twoSpin/2) );
  if( m_twoSpin % 2 == 1 ) return std::to_string( int(m_twoSpin/2) ) + "/2";
  return "";
}

ParticleProperties::ParticleProperties( const std::string& pdg_string )
{
  auto is_fundamental = []( unsigned pdgID)
  {
    return pdgID == 22 or pdgID == 24 or pdgID == 23 or (pdgID >= 11 and pdgID <=16 );
  };
  m_isValid = false;
  if ( pdg_string == "" || pdg_string.empty() ) return;
  if ( pdg_string[0] == '*' ) return;
  auto s = split( pdg_string, ',', false );
  if ( s.size() != 18 ) {
    DEBUG( "Invalid line : " << pdg_string );
    return;
  }
  std::transform(s.begin(), s.end(), s.begin(), trim);
  bool status    = 1;
  m_mass         = lexical_cast<double>( s[0], status );
  m_mErrPlus     = lexical_cast<double>( s[1], status );
  m_mErrMinus    = lexical_cast<double>( s[2], status );
  m_width        = lexical_cast<double>( s[3], status );
  m_wErrPlus     = lexical_cast<double>( s[4], status );
  m_wErrMinus    = lexical_cast<double>( s[5], status );
  m_pdgID        = lexical_cast<int>( s[12], status );
  m_Rexist       = lexical_cast<int>( s[14], status );
  m_Gparity      = chargeFromString( s[7], status );
  m_parity       = chargeFromString( s[9], status );
  m_Cparity      = chargeFromString( s[10], status );
  m_charge       = chargeFromString( s[13], status );
  m_isospin      = s[6];
  m_status       = s[15][0];
  m_name         = s[16];
  m_Aformat      = s[11].size() == 1 ? s[11][0] : ' ';
  m_chargeString = s[13];
  m_quarkContent = QuarkContent( s[17] );
  if( !is_fundamental(m_pdgID) && (s[17] == "??" or s[17] == "non-qQ" or  s[17] == "" ) )
  {
    WARNING("Quark content for: " << m_name << " not valid -> will be assumed to decay weakly"); 
  }

  bool spin_status = 1;
  if( s[8] == "?" ) m_twoSpin = -1;
  else if( s[8].find("/") != std::string::npos ){
    m_twoSpin = lexical_cast<int>( s[8].substr(0, s[8].find("/") ) , spin_status );
  }
  else m_twoSpin = 2 * lexical_cast<int>( s[8], spin_status );
  if( spin_status == 0 ){
    DEBUG("Spin of particle: " << name() << " could not be interpretted (J=" << s[8] << ")"  );
  }
  m_radius = 1.5 / GeV; 
  bool isCharm = ( abs(pdgID()) == 421 || 
                   abs(pdgID()) == 411 || 
                   abs(pdgID()) == 431 || 
                   abs(pdgID()) == 4122 );
  if(isCharm) m_radius = 5.0 / GeV; 
  m_isValid = true;
}

bool ParticleProperties::hasDistinctAnti() const { return !( m_Aformat == ' ' ); }

bool ParticleProperties::antiThis()
{
  if ( !hasDistinctAnti() ) return false;
  swapChars( m_chargeString, '+', '-');
  m_charge *= -1;
  if( isFermion() ) m_parity *= -1;
  m_quarkContent.antiThis();
  m_pdgID *= -1;
  return true;
}

ParticleProperties ParticleProperties::anti() const
{
  ParticleProperties PP( *this );
  PP.antiThis();
  return PP;
}

std::string ParticleProperties::name() const
{
  if( m_customName ) return m_name; 
  std::string fullName = m_name;
  if ( m_pdgID < 0 && m_Aformat == 'F' ) fullName += "bar";
  fullName += m_chargeString;
  return fullName;
}
double ParticleProperties::radius() const { return m_radius; }

double ParticleProperties::lifetime() const 
{ 
  return 6.582119514 / ( m_width * std::pow( 10, 13 ) );
}


bool ParticleProperties::isNonResonant() const
{
  return m_name.find("NonRes") != std::string::npos ; 
}

bool ParticleProperties::isFermion() const {
  return m_twoSpin % 2 == 1;
}

bool ParticleProperties::isPhoton() const 
{
  return m_pdgID == 22;
}

bool ParticleProperties::isNeutrino() const 
{
  return abs(m_pdgID) == 12 || abs(m_pdgID) == 14 || abs(m_pdgID) == 16;
}

bool ParticleProperties::isBoson() const {
  return ! isFermion();
}

std::ostream& AmpGen::operator<<( std::ostream& out, const ParticleProperties& pp )
{
  pp.print( out );
  return out;
}

std::string ParticleProperties::spinName() const {
  if( m_twoSpin % 2 == 0 ){
    if( m_twoSpin == 0 ) return "S";
    if( m_twoSpin == 2 ) return "V";
    if( m_twoSpin == 4 ) return "T";
    if( m_twoSpin == 6 ) return "H"; 
  }
  else {
    if( m_twoSpin == 1 ) return m_pdgID > 0 ? "f" : "F";
    if( m_twoSpin == 3 ) return m_pdgID > 0 ? "r" : "R";
    if( m_twoSpin == 5 ) return m_pdgID > 0 ? "e" : "E";
    if( m_twoSpin == 7 ) return m_pdgID > 0 ? "c" : "C";
  }
  WARNING("Spin name not implemented for " << m_twoSpin );
  return "?";
}

const ParticleProperties* ParticleProperties::get( const std::string& name, const bool& quiet )
{
  return ParticlePropertiesList::get( name, quiet );
}

bool ParticleProperties::operator==( const ParticleProperties& rhs ) const
{
  if ( pdgID() == 0 && rhs.pdgID() == 0 ) {
    return name() == rhs.name();
  }
  return pdgID() == rhs.pdgID();
}

bool ParticleProperties::operator<( const ParticleProperties& rhs ) const
{
  if ( pdgID() == 0 || rhs.pdgID() == 0 ) {
    return name() < rhs.name();
  }
  return pdgID() < rhs.pdgID();
}

bool ParticleProperties::operator>( const ParticleProperties& rhs ) const { return !( *this == rhs || *this < rhs ); }
bool ParticleProperties::operator<=( const ParticleProperties& rhs ) const { return ( *this < rhs || *this == rhs ); }
bool ParticleProperties::operator>=( const ParticleProperties& rhs ) const { return ( *this > rhs || *this == rhs ); }


void ParticleProperties::setProperty(const std::string& key, const std::string& value)
{
  DEBUG("Setting property: " << key << " " << value );
  bool status = true; 
  if( key == "mass" )  m_mass  = stod(value) / MeV;
  else if( key == "name" )  m_name = value; 
  else if( key == "width" ) m_width = stod(value) / MeV;
  else if( key == "spin"  ) 
  {
    bool spin_status = 1;
    if( value == "?" ) m_twoSpin = 0;
    else if( value.find("/") != std::string::npos ){
      m_twoSpin = lexical_cast<int>( value.substr(0, value.find("/") ) , spin_status );
    }
    else m_twoSpin = 2 * lexical_cast<int>( value, spin_status );
    if( spin_status == 0 ){
      ERROR("Spin of particle: " << name() << " could not be interpretted (J=" << value << ")"  );
    }
  }
  else if( key == "parity" ) m_parity = chargeFromString( value, status );
  else if( key == "charge" ) m_charge = chargeFromString( value, status );
  else if( key == "quarks" ) m_quarkContent = QuarkContent(value);
  else {
    ERROR("Unrecognised key: " << key );
  }
}
