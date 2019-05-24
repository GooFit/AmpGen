#include <stddef.h>
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "AmpGen/QuarkContent.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

std::array<char,6> QuarkState::gNames = {'d', 'u', 's', 'c', 'b', 't'};
std::map<char, int> QuarkState::gPositions;

QuarkState::QuarkState()
{
  for ( auto& quark : m_quarks ) quark = 0;
  initPositions();
}

QuarkState::QuarkState( const std::string& str ) : QuarkState() {
  for ( auto& c : str ) {
    auto lc = std::tolower(c);
    auto pos = gPositions.find( lc );
    if( pos == gPositions.end() ) continue; 
    m_quarks[ pos->second ] += ( islower(c) ? 1 : -1 );
  }
}

bool QuarkState::initPositions()
{
  if ( gPositions.empty() ) {
    for ( int i = 0; i < 6; i++ ) gPositions[gNames[i]] = i;
  }
  return true;
}

void QuarkState::antiThis()
{
  for ( auto& quark : m_quarks ) quark *= -1;
}

char QuarkState::nameFromPosition( int i ) const
{
  return ( i < 0 || i >= 6 ) ? 'X' : gNames[i];
}
int QuarkState::positionFromName( char c ) const
{
  std::map<char, int>::const_iterator it = gPositions.find( c );
  return gPositions.end() == it ? -9999 : it->second; 
}

void QuarkState::print( std::ostream& os ) const
{
  for ( unsigned int i = 0; i < m_quarks.size(); i++ ) {
    os << "(" << nameFromPosition( i ) << ":" << m_quarks[i] << ")";
    if ( i + 1 < m_quarks.size() ) os << " ";
  }
}

QuarkState& QuarkState::operator+=( const QuarkState& rhs )
{
  for ( unsigned int i = 0; i < 6; i++ ) m_quarks[i] += rhs[i];
  return *this;
}
QuarkState& QuarkState::operator-=( const QuarkState& rhs )
{
  for ( unsigned int i = 0; i < 6; i++ ) m_quarks[i] -= rhs[i];
  return *this;
}

QuarkState QuarkState::operator+( const QuarkState& rhs ) const
{
  QuarkState returnVal( *this );
  returnVal += rhs;
  return returnVal;
}
QuarkState QuarkState::operator-( const QuarkState& rhs ) const
{
  QuarkState returnVal( *this );
  returnVal -= rhs;
  return returnVal;
}
std::ostream& operator<<( std::ostream& st, const QuarkState& qc )
{
  qc.print( st );
  return st;
}
bool QuarkState::isVacuum() const
{
  for ( auto& quark : m_quarks )
    if ( quark != 0 ) return false;
  return true;
}

bool QuarkState::operator==( const QuarkState& rhs ) const
{
  for ( unsigned int i = 0; i < 6; ++i )
    if ( rhs[i] != m_quarks[i] ) return false;
  return true;
}

bool QuarkState::operator!=( const QuarkState& rhs ) const 
{
  return !( *this == rhs );
}

int QuarkState::operator[]( const size_t& index ) const { return m_quarks[index]; }

QuarkContent::QuarkContent() : m_quarks(1) {}

void QuarkContent::antiThis()
{
  for ( auto& qc : m_quarks ) qc.antiThis();
}

void QuarkContent::initFromString( const std::string& str )
{
  if ( str.find( "non-qQ" ) < str.size() ) {
    m_quarks.resize( 1 );
    return;
  }
  auto tokens = split( replaceAll( str, "sqrt", "" ), {'(', ')', '+', '-'} );
  m_quarks.clear();
  for ( auto& token : tokens ) {
    QuarkState qc( token );
    if ( !qc.isVacuum() ) m_quarks.emplace_back( qc );
  }
  if ( m_quarks.size() == 0 ) m_quarks.resize( 1 );
}

void QuarkContent::print( std::ostream& os ) const
{
  os << "[";
  for ( unsigned int i = 0; i < m_quarks.size(); i++ ) {
    m_quarks[i].print( os );
    if ( i != m_quarks.size() - 1 ) os << ", ";
  }
  os << "]";
}

QuarkContent& QuarkContent::operator+=( const QuarkContent& rhs )
{
  if ( size() * rhs.size() <= 1 ) {
    for ( unsigned int i = 0; i < size(); i++ ) m_quarks[i] += rhs[i];
   } else {
    QuarkContent oldThis( *this );
    m_quarks.clear();
    for ( unsigned int i = 0; i < oldThis.size(); i++ ) {
      for ( unsigned int j = 0; j < rhs.size(); j++ ) {
        m_quarks.push_back( oldThis[i] + rhs[j] );
      }
    }
  }
  return *this;
}
QuarkContent& QuarkContent::operator-=( const QuarkContent& rhs )
{
  if ( m_quarks.size() * rhs.quarks().size() <= 1 ) {
    for ( unsigned int i = 0; i < m_quarks.size(); i++ ) m_quarks[i] -= rhs[i];
  } else {
    QuarkContent oldThis( *this );
    m_quarks.clear();
    for ( unsigned int i = 0; i < oldThis.size(); i++ ) {
      for ( unsigned int j = 0; j < rhs.size(); j++ ) {
        m_quarks.push_back( oldThis[i] - rhs[j] );
      }
    }
  }
  return *this;
}

QuarkContent QuarkContent::operator+( const QuarkContent& rhs ) const
{
  QuarkContent returnVal( *this );
  returnVal += rhs;
  return returnVal;
}
QuarkContent QuarkContent::operator-( const QuarkContent& rhs ) const
{
  QuarkContent returnVal( *this );
  returnVal -= rhs;
  return returnVal;
}
std::ostream& AmpGen::operator<<( std::ostream& st, const QuarkContent& qc )
{
  qc.print( st );
  return st;
}

bool QuarkContent::operator==( const QuarkContent& rhs ) const
{
//  for ( auto& l : m_quarks ) {
//    if( std::any_of( rhs.quarks().begin(), rhs.quarks().end(), [&l](auto& r){ return l == r ; } ) ) 
//      return true; 
//  }
  for ( auto& l : m_quarks ) {
    for ( auto& r : rhs.quarks() ) {
      if ( l == r ) return true;
    }
  }
  return false;
}

bool QuarkContent::operator!=( const QuarkContent& rhs ) const 
{
  return !( *this == rhs );
}

QuarkState QuarkContent::operator[]( const size_t& index ) const 
{
  return m_quarks[index]; 
}

size_t QuarkContent::size() const { return m_quarks.size(); }

std::vector<QuarkState> QuarkContent::quarks() const 
{
  return m_quarks; 
}

