#ifndef AMPGEN_QUARKCONTENT_H
#define AMPGEN_QUARKCONTENT_H

#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace AmpGen
{
  class QuarkContent
  {
    static char _names[6];
    static char _NAMES[6];
    static std::map<char, int> _positions;
    std::array<int, 6> m_quarks;
    static bool initPositions();

  public:
    QuarkContent();
    QuarkContent( const std::string& str );

    void antiThis();
    bool initFromString( const std::string& str );

    char nameFromPosition( int i ) const;
    int positionFromName( char c ) const;
    bool isVacuum() const;

    void print( std::ostream& os = std::cout ) const;

    QuarkContent& operator+=( const QuarkContent& rhs );
    QuarkContent& operator-=( const QuarkContent& rhs );

    QuarkContent operator+( const QuarkContent& rhs ) const;
    QuarkContent operator-( const QuarkContent& rhs ) const;
    int operator[]( const size_t& index ) const;
    bool operator==( const QuarkContent& rhs ) const;
  };
} // namespace AmpGen
std::ostream& operator<<( std::ostream& st, const AmpGen::QuarkContent& qc );
#endif
//
