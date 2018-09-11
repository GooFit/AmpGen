#include <array>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "AmpGen/QuarkContent.h"

using namespace AmpGen;

char QuarkContent::_names[6] = {'d', 'u', 's', 'c', 'b', 't'};
char QuarkContent::_NAMES[6] = {'D', 'U', 'S', 'C', 'B', 'T'};
std::map<char, int> QuarkContent::_positions;

QuarkContent::QuarkContent()
{
  for ( auto& quark : m_quarks ) quark = 0;
  initPositions();
}

QuarkContent::QuarkContent( const std::string& str ) { initFromString( str ); }

bool QuarkContent::initPositions()
{
  if ( _positions.empty() ) {
    for ( int i = 0; i < 6; i++ ) _positions[_names[i]] = i;
    for ( int i = 0; i < 6; i++ ) _positions[_NAMES[i]] = i;
  }
  return true;
}

void QuarkContent::antiThis()
{
  for ( auto& quark : m_quarks ) quark *= -1;
}

bool QuarkContent::initFromString( const std::string& str )
{
  for ( auto& quark : m_quarks ) quark = 0;

  for ( unsigned int i = 0; i < str.size(); i++ ) {
    switch ( str[i] ) {
    case 'd':
      m_quarks[0]++;
      break;
    case 'D':
      m_quarks[0]--;
      break;
    case 'u':
      m_quarks[1]++;
      break;
    case 'U':
      m_quarks[1]--;
      break;
    case 's':
      m_quarks[2]++;
      break;
    case 'S':
      m_quarks[2]--;
      break;
    case 'c':
      m_quarks[3]++;
      break;
    case 'C':
      m_quarks[3]--;
      break;
    case 'b':
      m_quarks[4]++;
      break;
    case 'B':
      m_quarks[4]--;
      break;
    case 't':
      m_quarks[5]++;
      break;
    case 'T':
      m_quarks[5]--;
      break;
    }
  }
  return true;
}

char QuarkContent::nameFromPosition( int i ) const
{
  if ( i < 0 || i >= 6 )
    return 'X';
  else
    return _names[i];
}
int QuarkContent::positionFromName( char c ) const
{
  std::map<char, int>::const_iterator it = _positions.find( c );
  if ( _positions.end() == it )
    return -9999;
  else
    return it->second;
}

void QuarkContent::print( std::ostream& os ) const
{
  for ( unsigned int i = 0; i < m_quarks.size(); i++ ) {
    os << "(" << nameFromPosition( i ) << ":" << m_quarks[i] << ")";
    if ( i + 1 < m_quarks.size() ) os << " ";
  }
}

QuarkContent& QuarkContent::operator+=( const QuarkContent& rhs )
{
  for ( unsigned int i = 0; i < 6; i++ ) m_quarks[i] += rhs[i];
  return *this;
}
QuarkContent& QuarkContent::operator-=( const QuarkContent& rhs )
{
  for ( unsigned int i = 0; i < 6; i++ ) m_quarks[i] -= rhs[i];
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
std::ostream& operator<<( std::ostream& st, const QuarkContent& qc )
{
  qc.print( st );
  return st;
}
bool QuarkContent::isVacuum() const
{
  for ( auto& quark : m_quarks )
    if ( quark != 0 ) return false;
  return true;
}

bool QuarkContent::operator==( const QuarkContent& rhs ) const
{
  for ( unsigned int i = 0; i < 6; ++i )
    if ( rhs[i] != m_quarks[i] ) return false;
  return true;
}
int QuarkContent::operator[]( const size_t& index ) const { return m_quarks[index]; }
