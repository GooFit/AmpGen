#ifndef AMPGEN_QUARKCONTENT_H
#define AMPGEN_QUARKCONTENT_H

#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace AmpGen
{
  class QuarkState
  {
    private: 
      static std::array<char,6>  gNames; 
      static std::map<char, int> gPositions;
      static bool initPositions();
      std::array<int, 6> m_quarks;
    public:
      QuarkState();
      QuarkState( const std::string& str );
      void antiThis();
      void initFromString( const std::string& str );
      char nameFromPosition( int i ) const;
      int positionFromName( char c ) const;
      bool isVacuum()                const;
      void print( std::ostream& os = std::cout ) const;
      QuarkState& operator+=( const QuarkState& rhs );
      QuarkState& operator-=( const QuarkState& rhs );
      QuarkState  operator+ ( const QuarkState& rhs ) const;
      QuarkState  operator- ( const QuarkState& rhs ) const;
      bool        operator==( const QuarkState& rhs ) const;
      bool        operator!=( const QuarkState& rhs ) const;
      int         operator[]( const size_t& index ) const;
  };

  class QuarkContent
  {
    private: 
      std::vector<QuarkState> m_quarks;
    public:
      QuarkContent();
      QuarkContent( const std::string& str ){ initFromString(str); }
      void antiThis();
      void initFromString( const std::string& str );
      void print( std::ostream& os = std::cout ) const;
      size_t size() const; 
      bool compatible( const QuarkContent& other ) const;
      QuarkContent& operator+=( const QuarkContent& rhs );
      QuarkContent& operator-=( const QuarkContent& rhs );
      QuarkContent  operator+ ( const QuarkContent& rhs ) const;
      QuarkContent  operator- ( const QuarkContent& rhs ) const;
      bool          operator==( const QuarkContent& rhs ) const;
      bool          operator!=( const QuarkContent& rhs ) const;
      QuarkState    operator[]( const size_t& index) const; 
      std::vector<QuarkState> quarks() const; 
  };
  std::ostream& operator<<( std::ostream& st, const QuarkState& qc );
  std::ostream& operator<<( std::ostream& st, const QuarkContent& qc );
} // namespace AmpGen
#endif
