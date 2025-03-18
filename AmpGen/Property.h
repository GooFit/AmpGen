#ifndef AMPGEN_PROPERTY_H
#define AMPGEN_PROPERTY_H
#include <stddef.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <map>
#include <cstring>

#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen {
  template <typename value_t> 
    class Property {
      protected:
        std::string m_name; 
        std::string m_helpString; 
        value_t m_value; 
      public: 
        Property( const std::string& name, const value_t& def=value_t(), const std::string& helpString="" ) : 
          m_name(name),
          m_helpString(helpString),
          m_value(def) {
            setFromOptionsParser();
            if ( OptionsParser::printHelp() ) help(def);
            DEBUG( *this );
          }

        template <class G> bool operator==(const G& other) const { return m_value == other; }
        template <class G> bool operator!=(const G& other) const { return m_value != other; }
        operator value_t() const { return m_value; }
        operator value_t()       { return m_value; }
        const std::string& name() const { return m_name ; } 
        template <typename T> friend std::ostream& operator<<( std::ostream& os, const Property<T>& np );
      private: 
        void help(const value_t& def){
          std::string type = type_string<value_t>();
          if( type == "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >" ) type = "string";
          std::cout << " " << bold_on << std::left << std::setw(27) << m_name << bold_off << std::setw(20) << "[" + type + "]" ;
          auto tokens = split( m_helpString, '\n' );
          if( tokens.size() == 0 ) std::cout << std::endl; 
          for( size_t i = 0 ; i < tokens.size(); ++i){
            if( i == 0 ){
              std::cout << tokens[i];
              if( def != value_t() ) std::cout << " (default = " << def << ")";
              std::cout << std::endl;  
            }
            else std::cout << std::string(48,' ') << tokens[i] << std::endl; 
          }
        }
        bool setFromOptionsParser(){
          auto parser = OptionsParser::getMe();
          auto line = parser->find( m_name );
          if( line == parser->end() ) return false ; 
          const std::vector<std::string>& vsl = line->second;
          if ( vsl.size() < 2 ) return false; // first element is parameter name
          bool status = true;
          if constexpr( isVector<value_t>::value ){
            m_value.resize( vsl.size() - 1 ); 
            for ( unsigned int i = 1; i < vsl.size(); i++ ) {
              m_value[i-1] = lexical_cast<value_t::value_type>( vsl[i], status );
              if ( status == false ) {
                ERROR( "Failed to parse token: " << vsl[i] << " for parameter: " << m_name );
                return false; 
              }
            }
          }
          else {
            if( vsl.size() != 2 ){
              ERROR("Constructing scalar quantity, only one argument expected, but " << vsl.size() -1 << " found"); 
              return false; 
            }
            m_value = lexical_cast<value_t>( vsl[1], status );
            if ( status == false ) {
              ERROR( "Failed to parse token: " << vsl[1] << " for parameter: " << m_name );
              return false; 
            }
          }
          return true;
        }
    };
    template <typename T> std::ostream& operator<<( std::ostream& os, const Property<T>& np );
    template <typename ...T>  std::string helpStringOptions(const std::string& header, const T&... args);
};

template <typename T> std::ostream& AmpGen::operator<<( std::ostream& os, const AmpGen::Property<T>& np ) {
  os << np.name() << " = ";
  if constexpr( ! isVector<T>::value ) os << np.m_value; 
  else {
    for( unsigned i{0} ; i != np.m_value.size() ; ++i ) os << np.m_value[i] << " ";
  }
  return os;
}

template <typename ...T> std::string AmpGen::helpStringOptions(const std::string& header, const T&... args )
{
  std::stringstream rt;
  rt << header;
  for_each( std::make_tuple(args...), [&rt](const auto& f) mutable {
    rt << "\n\033[3m "  << f.first << "\033[0m: " << f.second; 
  });
  return rt.str();
}
#endif
