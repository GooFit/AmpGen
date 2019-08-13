#include "AmpGen/MsgService.h"

#define declare_enum(name, ...) enum name {__VA_ARGS__};   \
template <class T = name> T parse(const std::string& word){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::parse<T>(word, args); } \
template <class T = name> std::string to_string( const T& enumItem ){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::to_string<T>(enumItem, args) ; }

namespace AmpGen {
  namespace detail {
    template <class T> T parse(const std::string& word, const char* args)
    {
      char* p;                                                 
      auto number = strtoul( word.c_str(), &p, 10 );           
      if( *p == 0 ) return T(number);
      size_t counter = 0;
      size_t begin = 0;
      size_t end   = 0;
      auto compare = [](const char* word, const char* otherWord, const size_t& nChar)
      {
        for( size_t x = 0; x != nChar ; ++x) if( word[x] != otherWord[x] ) return false;
        return true;
      };
      for( ; args[begin] != '\0' ; begin++ )
      {
        while( args[begin] == ' ' ) begin++;
        for( end=begin; args[end] != '\0'; end++ ) if( args[end] == ',' ) break;
        if( compare( word.c_str(), args + begin , end-begin ) ) break; 
        begin = end;
        counter++; 
      }
      if( args[begin] == '\0' ) return T(counter-1);
      return T(counter);                                       
    }
    template <class T> std::string to_string(const T& enumItem, const char* args)
    {
      size_t counter = 0;
      size_t sBegin  = 0;
      size_t sLength = 0;
      for( ; args[sBegin] != '\0' && counter != enumItem; sBegin++ )
      {
        if( args[sBegin] == ',' ) counter++;
      }
      while( args[sBegin] == ' ') sBegin++;
      for(; args[sLength + sBegin] != '\0' ; ++sLength ) if( args[sBegin + sLength] == ',') break; 
      return std::string(args).substr(sBegin, sLength);
    }
  }
}
