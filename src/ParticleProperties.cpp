// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:04 GMT
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"

#include <algorithm>
#include <stdlib.h>
#include <vector>

#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

double ParticleProperties::_defaultRadius      = 1.5 / GeV; // 1.5/GeV;
double ParticleProperties::_defaultCharmRadius = 5.0 / GeV;

void ParticleProperties::print( std::ostream& out ) const
{
  out << "Mass  " << mass() << " +" << mErrPlus() << " -" << mErrMinus() << "\nWidth " << width() << " +" << wErrPlus()
      << " -" << wErrMinus() << "\n I=" << I() << ", G=" << G() << "\n J=" << J() << ", C=" << C() << ", P=" << P()
      << "\n Q = " << charge() << "\n pdgID " << pdgID() << "\n name " << name() << "\n quark content: " << quarks()
      << "\n net-quark-content " << netQuarkContent() << "\n is its own antiparticle? "
      << ( isItsOwnAnti() ? "yes" : "no" ) << "\n radius " << radius() * GeV << " /GeV"
      << "\n";
}

int ParticleProperties::chargeFromString( const std::string& ch, bool& status ) const
{

  if ( ch == "+" ) return 1;
  if ( ch == "-" ) return -1;
  if ( ch == " " ) return 0;
  if ( ch == "0" ) return 0;
  if ( ch == "++" ) return 2;
  if ( ch == "--" ) return -2;
  status = 0;
  return 0;
}

ParticleProperties::ParticleProperties( const std::string& pdg_string ) : m_netQuarkContent()
{
  m_isValid = false;
  if ( pdg_string == "" || pdg_string.empty() ) return;
  if ( pdg_string[0] == '*' ) return;
  auto s = split( pdg_string, ',', false );
  if ( s.size() != 18 ) {
    DEBUG( "Invalid line : " << pdg_string );
    return;
  }
  for ( auto& st : s ) st = trim( st );
  bool status             = 1;
  m_mass                  = lexical_cast<double>( s[0], status );
  m_mErrPlus              = lexical_cast<double>( s[1], status );
  m_mErrMinus             = lexical_cast<double>( s[2], status );
  m_width                 = lexical_cast<double>( s[3], status );
  m_wErrPlus              = lexical_cast<double>( s[4], status );
  m_wErrMinus             = lexical_cast<double>( s[5], status );
  m_pdgID                = lexical_cast<int>( s[12], status );
  m_Rexist                = lexical_cast<int>( s[14], status );

  m_Gparity      = chargeFromString( s[7], status );
  m_Parity       = chargeFromString( s[9], status );
  m_Cparity      = chargeFromString( s[10], status );
  m_charge       = chargeFromString( s[13], status );
  m_Isospin      = s[6];
  m_JtotalSpin   = s[8];
  m_status       = s[15][0];
  m_name         = s[16];
  m_quarks       = s[17];
  m_Aformat      = s[11][0];
  m_chargeString = s[13];
  m_netQuarkContent.initFromString( m_quarks );
  bool spin_status = 1;
  if( m_JtotalSpin == "?" ) m_twoSpin = 0;
    else if( m_JtotalSpin.find("/") != std::string::npos ){
    m_twoSpin = lexical_cast<int>( m_JtotalSpin.substr(0, m_JtotalSpin.find("/") ) , spin_status );
  }
  else m_twoSpin = 2 * lexical_cast<int>( m_JtotalSpin, spin_status );
  if( spin_status == 0 ){
   DEBUG("Spin of particle: " << name() << " could not be interpretted (J=" << m_JtotalSpin << ")"  );
  }
  setRadius();

  m_isValid = true;
}

void ParticleProperties::setRadius()
{
  // set radius (not part of mass_width.csv):
  bool isCharm = ( abs( pdgID() ) == 421 || abs( pdgID() ) == 411 || abs( pdgID() ) == 431 );
  m_Radius     = isCharm ? _defaultCharmRadius : _defaultRadius;
}

void ParticleProperties::antiQuarks()
{
  if ( m_quarks.empty() ) return;

  replace( m_quarks.begin(), m_quarks.end(), 'U', 'h' );
  replace( m_quarks.begin(), m_quarks.end(), 'D', 'l' );
  replace( m_quarks.begin(), m_quarks.end(), 'C', 'i' );
  replace( m_quarks.begin(), m_quarks.end(), 'S', 'm' );
  replace( m_quarks.begin(), m_quarks.end(), 'T', 'j' );
  replace( m_quarks.begin(), m_quarks.end(), 'B', 'n' );

  replace( m_quarks.begin(), m_quarks.end(), 'u', 'U' );
  replace( m_quarks.begin(), m_quarks.end(), 'd', 'D' );
  replace( m_quarks.begin(), m_quarks.end(), 'c', 'C' );
  replace( m_quarks.begin(), m_quarks.end(), 's', 'S' );
  replace( m_quarks.begin(), m_quarks.end(), 't', 'T' );
  replace( m_quarks.begin(), m_quarks.end(), 'b', 'B' );

  replace( m_quarks.begin(), m_quarks.end(), 'h', 'u' );
  replace( m_quarks.begin(), m_quarks.end(), 'l', 'd' );
  replace( m_quarks.begin(), m_quarks.end(), 'i', 'c' );
  replace( m_quarks.begin(), m_quarks.end(), 'm', 's' );
  replace( m_quarks.begin(), m_quarks.end(), 'j', 't' );
  replace( m_quarks.begin(), m_quarks.end(), 'n', 'b' );

  unsigned int pos = m_quarks.find( "SqrT" );
  if ( pos < m_quarks.size() ) {
    m_quarks.replace( pos, 4, "sqrt" );
  }
}

void ParticleProperties::antiQuarkContent() { m_netQuarkContent.antiThis(); }

void ParticleProperties::antiCharge()
{
  replace( m_chargeString.begin(), m_chargeString.end(), '+', 'f' );
  replace( m_chargeString.begin(), m_chargeString.end(), '-', '+' );
  replace( m_chargeString.begin(), m_chargeString.end(), 'f', '-' );
  m_charge *= -1;
}
bool ParticleProperties::hasDistinctAnti() const { return !( m_Aformat == ' ' ); }

bool ParticleProperties::barred() const { return m_Aformat == 'F'; }

bool ParticleProperties::antiThis()
{
  if ( !hasDistinctAnti() ) return false;
  antiCharge();
  antiQuarks();
  antiQuarkContent();
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
  std::string fullName = m_name;
  if ( m_pdgID < 0 && m_Aformat == 'F' ) fullName += "bar";
  fullName += m_chargeString;
  return fullName;
}
double ParticleProperties::radius() const { return m_Radius; }

bool ParticleProperties::isNonResonant() const
{
  return m_name.find("NonRes") != std::string::npos ; 
}

bool ParticleProperties::isFermion() const {
  return m_twoSpin % 2 == 1;
}

bool ParticleProperties::isBoson() const {
  return ! isFermion();
}

std::ostream& operator<<( std::ostream& out, const ParticleProperties& pp )
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
   if( m_twoSpin == 1 ) return "f";
   if( m_twoSpin == 3 ) return "r";
   if( m_twoSpin == 5 ) return "e";
   if( m_twoSpin == 7 ) return "c";
  }
  WARNING("Spin name not implemented for " << m_JtotalSpin );
  return "?";
}

const ParticleProperties* ParticleProperties::get( const std::string& name, const bool& quiet ){
  return ParticlePropertiesList::get( name, quiet );
}

