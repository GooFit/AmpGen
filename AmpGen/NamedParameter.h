#ifndef AMPGEN_NAMEDPARAMETER_H
#define AMPGEN_NAMEDPARAMETER_H
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:56 GMT

#include <stddef.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <map>

#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen
{
  /** @class NamedParameter
      @brief A parameter with value specified by the user at runtime, either in an options file or via the command line

      Stores a vector of values for a parameter.
      @tparam T the type of this named parameter, i.e. strings, or numbers or bools etc. 

    */
  template <class T>
  class NamedParameter 
  {
  protected:
    std::string    m_name;         /// < Name of this parameter
    std::string    m_helpString;   /// < The helper string for this parameter, printed if the flag --help is used. 
    std::vector<T> m_valueArray;   /// < The value (array) of this parameter. 

    bool setFromOptionsParser()
    {
      auto parser = OptionsParser::getMe();
      auto line = parser->find( m_name );
      if( line == parser->end() ) return false ; 
      const std::vector<std::string>& vsl = line->second;
      if ( vsl.size() < 2 ) return false; // first element is parameter name
      m_valueArray.clear();
      m_valueArray.resize( vsl.size() - 1 );
      for ( unsigned int i = 1; i < vsl.size(); i++ ) {
        bool status = true;
        T tmpVal    = lexical_cast<T>( vsl[i], status );
        if ( status == false ) {
          ERROR( "Failed to parse token: " << vsl[i] << " for parameter: " << m_name );
          continue;
        }
        setVal( tmpVal, i - 1 );
      }
      return true;
    }

  public:
    NamedParameter( const std::string& name, const T& def=T(), const std::string& helpString="" ) : 
      m_name(name),
      m_helpString(helpString)
    {
      setVal( def );
      setFromOptionsParser();
      if ( OptionsParser::printHelp() ) help(def);
      DEBUG( *this );
    }
    NamedParameter(const std::string& name, const std::vector<T>& defVec, const std::string& helpString="")
        : m_name(name),
          m_helpString(helpString)
    {
      setVal( defVec );
      setFromOptionsParser();
      if ( OptionsParser::printHelp() ) help( defVec.size() > 0 ? defVec[0] : T() );
    }
    
    void help(const T& def){
      std::map< std::string, std::string > aliases;  
      std::string type = typeof<T>();
      if( type == "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >" ) type = "string";
      std::cout << " " << bold_on << std::left << std::setw(27) << m_name << bold_off 
       << std::setw(20) << "[" + type + "]" ;
      auto tokens = split( m_helpString, '\n' );
      if( tokens.size() == 0 ) std::cout << std::endl; 
      for( size_t i = 0 ; i < tokens.size(); ++i){
        if( i == 0 ){
          std::cout << tokens[i];
          if( def != T() ) std::cout << " (default = " << def << ")";
          std::cout << std::endl;  
        }
        else std::cout << std::string(48,' ') << tokens[i] << std::endl; 
      }
    }

    size_t size() const { return m_valueArray.size(); }

    const T getVal( int i = 0 ) const 
    {
      if ( i < 0 || i >= (int)m_valueArray.size() ) {
        ERROR( "Parameter name \"" << name() << "\". Index: " << i << " out of range [ " << 0 << " , "
                                   << m_valueArray.size() );
        throw std::runtime_error( "array index out of bounds" );
      }
      return m_valueArray[i];
    }

    operator T() const { return getVal(); }
    operator T()       { return getVal(); }
    const std::vector<T>& getVector() const { return m_valueArray; }

    void setVal( const T& val, int i = 0 )
    {
      if ( i < 0 ) return;
      if ( i >= (int)m_valueArray.size() ) {
        m_valueArray.resize( i + 1 );
      }
      m_valueArray[i]  = val;
    }
    void setVal( const std::vector<T>& valList )
    {
      m_valueArray     = valList;
    }

    operator std::vector<T>() const { return getVector(); }
    NamedParameter<T>& operator=( const T& d )
    {
      setVal( d, 0 );
      return *this;
    }
    NamedParameter<T>& operator=( const std::vector<T>& v )
    {
      setVal( v );
      return *this;
    }
    const std::string& name() const { return m_name ; } 

    static std::vector<T> getVectorArgument( const std::string& name, const T& default_value )
    {
      std::vector<T> return_container;
      unsigned int x = 0;
      T obj          = default_value;
      do {
        obj = AmpGen::NamedParameter<T>( name + std::to_string( x++ ), default_value );
        if ( obj != T() ) return_container.push_back( obj );
      } while ( obj != default_value );
      return return_container;
    }
  };
  template <class T> 
  std::ostream& operator<<( std::ostream& os, const NamedParameter<T>& np );
  
  std::string optionalHelpString(const std::string& header, const std::vector<std::pair<std::string, std::string>>& args);
}


template <typename T>
std::ostream& AmpGen::operator<<( std::ostream& os, const AmpGen::NamedParameter<T>& np )
{
  os  << np.name() ;
  for ( size_t i = 0; i < np.size(); i++ ) {
    if( i == 0 ) os << " = ";
    os << np.getVal( i );
    if ( i != np.size() ) os << " ";
  }
  return os;
}

#endif
//
