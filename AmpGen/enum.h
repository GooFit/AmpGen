#ifndef AMPGEN_ENUM_H
#define AMPGEN_ENUM_H 1
#include "AmpGen/MsgService.h"

#define declare_enum(name, ...)  \
enum class name {__VA_ARGS__};                                     \
template <> name parse(const std::string& word);                   \
template <> std::string to_string( const name& enumItem );         \
std::ostream& operator<<( std::ostream& os, const name& np); 

#define complete_enum(name, ...)                 \
template <> name parse(const std::string& word){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::parse<name>(word, args); } \
template <> std::string to_string( const name& enumItem ){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::to_string<name>(enumItem, args) ; } \
template <> name lexical_cast(const std::string& word, bool& /*status*/){ return parse<name>(word); } \
std::ostream& operator<<(std::ostream& os, const name& np){ return os << to_string<name>(np);}

#define make_enum(name, ...)                                                                                                                   \
enum class name {__VA_ARGS__};                                                                                                                 \
template <> name parse(const std::string& word){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::parse<name>(word, args); } \
template <> std::string to_string( const name& enumItem ){ constexpr auto args = #__VA_ARGS__; return AmpGen::detail::to_string<name>(enumItem, args) ; } \
template <> name lexical_cast(const std::string& word, bool& /*status*/){ return parse<name>(word); } \
std::ostream& operator<<(std::ostream& os, const name& np){ return os << to_string<name>(np);}

namespace AmpGen {
  template <class T> T parse( const std::string& word ){ return T(); }
  template <class T> std::string to_string( const T& enumItem ){ return ""; }

  namespace detail {
    template <class T> T parse(const std::string& word, const char* args)
    {
      char* p;                                                 
      auto number = strtoul( word.c_str(), &p, 10 );           
      if( *p == 0 ) return T(number);
      unsigned counter = 0;
      unsigned begin = 0;
      unsigned end   = 0;
      auto compare = [](const char* word, const char* otherWord, const unsigned& nChar)
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
      unsigned counter = 0;
      unsigned sBegin  = 0;
      unsigned sLength = 0;
      for( ; args[sBegin] != '\0' && counter != unsigned(enumItem); sBegin++ )
      {
        if( args[sBegin] == ',' ) counter++;
      }
      while( args[sBegin] == ' ') sBegin++;
      for(; args[sLength + sBegin] != '\0' ; ++sLength ) if( args[sBegin + sLength] == ',') break; 
      return std::string(args).substr(sBegin, sLength);
    }
  }
}
#endif
