#include "AmpGen/MultiQuarkContent.h"

#include <ext/alloc_traits.h>
#include <memory>

#include "AmpGen/Utilities.h"

using namespace AmpGen;

MultiQuarkContent::MultiQuarkContent() : std::vector<QuarkContent>( 1 ) {}

MultiQuarkContent::MultiQuarkContent( const MultiQuarkContent& other )

    = default;

void MultiQuarkContent::antiThis()
{
  for ( auto& qc : *this ) qc.antiThis();
}

bool MultiQuarkContent::initFromString( const std::string& str )
{

  if ( str.find( "non-qQ" ) < str.size() ) {
    std::vector<QuarkContent>::resize( 1 );
    return true;
  }
  auto tokens = split( replaceAll( str, "sqrt", "" ), {'(', ')', '+', '-'} );
  std::vector<QuarkContent>::clear();

  for ( auto& token : tokens ) {
    QuarkContent qc( token );
    if ( !qc.isVacuum() ) std::vector<QuarkContent>::emplace_back( qc );
  }
  if ( size() == 0 ) std::vector<QuarkContent>::resize( 1 );

  return true;
}

void MultiQuarkContent::print( std::ostream& os ) const
{
  os << "[";
  for ( unsigned int i = 0; i < this->size(); i++ ) {
    ( *this )[i].print( os );
    if ( i != size() - 1 ) os << ", ";
  }
  os << "]";
}

MultiQuarkContent& MultiQuarkContent::operator+=( const MultiQuarkContent& rhs )
{
  if ( this->size() * rhs.size() <= 1 ) {
    for ( unsigned int i = 0; i < this->size(); i++ ) ( *this )[i] += rhs[i];
  } else {
    MultiQuarkContent oldThis( *this );
    this->clear();
    for ( unsigned int i = 0; i < oldThis.size(); i++ ) {
      for ( unsigned int j = 0; j < rhs.size(); j++ ) {
        this->push_back( oldThis[i] + rhs[j] );
      }
    }
  }
  return *this;
}
MultiQuarkContent& MultiQuarkContent::operator-=( const MultiQuarkContent& rhs )
{
  if ( this->size() * rhs.size() <= 1 ) {
    for ( unsigned int i = 0; i < this->size(); i++ ) ( *this )[i] -= rhs[i];
  } else {
    MultiQuarkContent oldThis( *this );
    this->clear();
    for ( unsigned int i = 0; i < oldThis.size(); i++ ) {
      for ( unsigned int j = 0; j < rhs.size(); j++ ) {
        this->push_back( oldThis[i] - rhs[j] );
      }
    }
  }
  return *this;
}

MultiQuarkContent MultiQuarkContent::operator+( const MultiQuarkContent& rhs ) const
{
  MultiQuarkContent returnVal( *this );
  returnVal += rhs;
  return returnVal;
}
MultiQuarkContent MultiQuarkContent::operator-( const MultiQuarkContent& rhs ) const
{
  MultiQuarkContent returnVal( *this );
  returnVal -= rhs;
  return returnVal;
}
std::ostream& AmpGen::operator<<( std::ostream& st, const MultiQuarkContent& qc )
{
  qc.print( st );
  return st;
}
bool MultiQuarkContent::operator==( const MultiQuarkContent& rhs ) const
{
  for ( auto& l : *this ) {
    for ( auto& r : rhs ) {
      if ( l == r ) return true;
    }
  }
  return false;
}
