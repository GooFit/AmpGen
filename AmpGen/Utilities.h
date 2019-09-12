#ifndef AMPGEN_UTILITIES_H
#define AMPGEN_UTILITIES_H

#include <cxxabi.h>
#include <stddef.h>
#include <algorithm>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include <future>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "AmpGen/MsgService.h"
#include "AmpGen/MetaUtils.h"
namespace AmpGen {
  template <class T>
    bool isIn( const std::vector<T>& container, const T& obj )
    {
      for ( auto& it : container )
        if ( obj == it ) return true;
      return false;
    }

  template <class T, class B, class F>
    bool isIn( const std::vector<T>& container, const B& obj, F f )
    {
      for ( auto& it : container )
        if ( f( it, obj ) ) return true;
      return false;
    }

  template <class T>
    std::string vectorToString( const std::vector<T>& obj, const std::string& delim = "" )
    {
      std::stringstream ss;
      if( obj.size() == 0 ) return "";
      for ( unsigned int i = 0 ; i < obj.size()-1; ++i ) 
        ss << obj[i] << delim;
      ss << obj[obj.size()-1];
      return ss.str();
    }

  template <class T, class F>
    std::string vectorToString( const std::vector<T>& obj, const std::string& delim="", const F& functor =[](const T& f){ return f ; }  )
    {
      std::stringstream ss;
      if( obj.size() == 0 ) return "";
      for ( unsigned int i = 0 ; i < obj.size()-1; ++i ) 
        ss << functor(obj[i]) << delim;
      ss << functor(obj[obj.size()-1]);
      return ss.str();
    }

  template <class T> std::vector<std::vector<T>> nCr( const T& n, const T& r )
    {
      std::vector<bool> mask( n );
      std::vector<std::vector<T>> combinations;
      std::fill( mask.begin() + r, mask.end(), true );
      do {
        std::vector<T> perm;
        for ( T i = 0; i < n; ++i ) if ( !mask[i] ) perm.push_back(i);
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

  unsigned int FNV1a_hash( const std::string& toHash );

  std::vector<std::string> getItems( const std::string& tree, const std::vector<std::string>& brackets = {"{", "}"},
      const std::string& seperator = "," );

  void swapChars(std::string& arg, const char a, const char b );

  unsigned int editDistance( const std::string& s1, const std::string& s2 );

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

  template <class ...ARGS>
    std::string mysprintf(const std::string& format,
        ARGS&&... args){
      auto size = std::snprintf(nullptr, 0, format.c_str(), std::forward<ARGS>(args)...);
      std::string output(size+1,'\0');
      std::sprintf(&output[0],format.c_str(), std::forward<ARGS>(args)...);
      return output.substr(0, output.size()-1);
    }


  template <> double       lexical_cast( const std::string& word, bool& status );
  template <> unsigned int lexical_cast( const std::string& word, bool& status );
  template <> std::string  lexical_cast( const std::string& word, bool& status );
  template <> float        lexical_cast( const std::string& word, bool& status );
  template <> bool         lexical_cast( const std::string& word, bool& status );
  template <> int          lexical_cast( const std::string& word, bool& status );
  template <> unsigned long int lexical_cast( const std::string& word, bool& status );
  template <> long int          lexical_cast( const std::string& word, bool& status );

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

  // calculate number of swaps required to reorder set, where swaps can have different exchange parities. 
  std::pair<size_t,int> minSwaps(const std::vector<size_t>& indices, const std::vector<int>& exchangeParities); 

  template<class iterator, 
    class comparator>
      void parallel_sort(iterator begin, 
          iterator end, 
          const comparator& comp,
          const size_t& grainsize)
      {
        const size_t len   = end - begin; 
        if(len < grainsize) std::sort(begin, end, comp);
        else
        {
          const auto  middle = begin + len/2;
          auto future = std::async(parallel_sort<iterator,comparator>, begin, middle, comp, grainsize);
          parallel_sort(middle, end, comp, grainsize);
          future.wait();
          std::inplace_merge(begin, middle, end, comp);
        }
      }
  
  template<class iterator, 
    class initial_value,
    class functor>
      initial_value parallel_accumulate(iterator begin, 
          iterator end, 
          const initial_value& init,
          const functor& f)
      {
        auto total = init; 
        auto size  = end-begin;
        #ifdef _OPENMP
        #pragma omp parallel for reduction( +: total )
        #endif
        for( int it = 0; it < size; ++it ){
          total += f(*(begin+it));
        }
        return total;
      }

  template<class iterator>
    void parallel_sort(iterator begin, 
        iterator end, 
        const size_t& grainsize)
    {
      typedef typename std::iterator_traits<iterator>::value_type value_type;
      parallel_sort(begin,end,std::less<value_type>(), grainsize);
    }

  bool stringMatchesWildcard( const std::string& input, const std::string& wildcard_string,
      const char wildcard_character = '*' );

  bool isDir( const std::string& fname );
  bool fileExists( const std::string& name );

  std::vector<std::string> getListOfFiles(const std::string& directory, const std::string& patternString = "");

  void printSplash();
  void printReleaseNotes(const std::string& fname);

  std::string ltrim( std::string s );
  std::string rtrim( std::string s );
  std::string trim( std::string s );
  std::string expandGlobals( std::string path );

  std::ostream& bold_on( std::ostream& );
  std::ostream& bold_off( std::ostream& );
  std::ostream& italic_on( std::ostream& );
  std::ostream& italic_off( std::ostream& );
}
#endif
