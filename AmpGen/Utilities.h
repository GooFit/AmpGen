#ifndef AMPGEN_UTILITIES_H
#define AMPGEN_UTILITIES_H

#include <algorithm>
#include <cstdint>
#include <cxxabi.h>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "AmpGen/MsgService.h"
#include "AmpGen/MetaUtils.h"
/*

   Generic utility functions for AmpGen library
   nCr( n , r ) - Gets all unique (1..N) choose r elements
   vectorToString<T> - concatenates vector into a string
   findAll(input, ch) - Finds all positions of ch input string input
   vecFindAll ( input, strings ) - Finds positions of elements in strings in input,
   puts them into map
   getItems - Extract vector of head / branches in a persistified decay tree

*/

template <class T>
static bool isIn( const std::vector<T>& container, const T& obj )
{
  for ( auto& it : container )
    if ( obj == it ) return true;
  return false;
}

template <class T, class B, class F>
static bool isIn( const std::vector<T>& container, const B& obj, F f )
{
  for ( auto& it : container )
    if ( f( it, obj ) ) return true;
  return false;
}

template <class T>
std::string vectorToString( const std::vector<T>& obj, const std::string& delim = "" )
{
  std::string returnValue;
  std::stringstream ss;
  if( obj.size() == 0 ) return "";
  for ( unsigned int i = 0 ; i < obj.size()-1; ++i ) 
    ss << obj[i] << delim;
  ss << obj[obj.size()-1];
  return ss.str();
}



template <class T> 
std::vector<std::vector<T>> nCr( const T& n, const T& r )
{

  std::vector<bool> mask( n );
  std::vector<std::vector<T>> combinations;

  std::fill( mask.begin() + r, mask.end(), true );
  do {
    std::vector<T> perm;
    for ( T i = 0; i < n; ++i ) {
      if ( !mask[i] ) perm.push_back( i );
    }
    combinations.push_back( perm );
  } while ( std::next_permutation( mask.begin(), mask.end() ) );
  return combinations;
}
std::vector<std::string> vectorFromFile( const std::string& filename, const char ignoreLinesThatBeginWith = '#' );

std::vector<std::string> split( const std::string& s, char delim, bool ignoreWhitespace = true );
std::vector<std::string> split( const std::string& s, const std::vector<char>& delims );

std::vector<size_t> findAll( const std::string& input, const std::string& ch );

std::map<size_t, std::string> vecFindAll( const std::string& input, const std::vector<std::string>& vCh );

size_t find_next_of( const std::string& input, const std::vector<std::string>& patterns, const size_t& begin = 0 );

std::string replaceAll( const std::string& input, const std::string& toReplace, const std::string& replaceWith );

std::string convertTeXtoROOT( std::string input );

unsigned int FNV1a_hash( const std::string& toHash );

std::vector<std::string> getItems( const std::string& tree, const std::vector<std::string>& brackets = {"{", "}"},
                                   const std::string& seperator = "," );

unsigned int edit_distance( const std::string& s1, const std::string& s2 );

std::string round( const double& number, const unsigned int& nsf );

std::string numberWithError( const double& number, const double& error, const unsigned int& nDigits );

template <class RETURN_TYPE>
RETURN_TYPE lexical_cast( const std::string& word, bool& status )
{
  WARNING( "Only use specialised versions of this template (word = " << word << ", type = " << AmpGen::typeof<RETURN_TYPE>()
                                                                     << ")  " );
  status = 0;
  return RETURN_TYPE();
}

template <>
double lexical_cast( const std::string& word, bool& status );
template <>
unsigned int lexical_cast( const std::string& word, bool& status );
template <>
std::string lexical_cast( const std::string& word, bool& status );
template <>
float lexical_cast( const std::string& word, bool& status );
template <>
bool lexical_cast( const std::string& word, bool& status );
template <>
int lexical_cast( const std::string& word, bool& status );
template <>
uint64_t lexical_cast( const std::string& word, bool& status );
template <>
int64_t lexical_cast( const std::string& word, bool& status );

template <class FCN>
void processFile( const std::string& filename, FCN&& toDo, const char ignoreLinesThatBeginWith = '#' )
{
  std::string tmp;
  std::ifstream inFile( filename.c_str() );
  while ( inFile.good() ) {
    std::getline( inFile, tmp );
    if ( tmp.size() == 0 || tmp[0] == ignoreLinesThatBeginWith ) continue;
    toDo( tmp );
  }
  inFile.close();
}

bool stringMatchesWildcard( const std::string& input, const std::string& wildcard_string,
                            const char wildcard_character = '*' );

std::vector<std::string> getListOfFiles( const std::string& directory, std::string patternString = "" );
void printSplash();

// trim from end
std::string ltrim( std::string s );
std::string rtrim( std::string s );
// trim from both ends
std::string trim( std::string s );

bool file_exists( const std::string& name );

std::string expandGlobals( std::string path );

std::ostream& bold_on( std::ostream& );
std::ostream& bold_off( std::ostream& );
std::ostream& italic_on( std::ostream& );
std::ostream& italic_off( std::ostream& );
#endif
