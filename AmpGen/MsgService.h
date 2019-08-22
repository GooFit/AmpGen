#ifndef AMPGEN_MSGSERVICE_H
#define AMPGEN_MSGSERVICE_H

/** @defgroup msgService Messaging and logging
   MsgService Header
   Defines coloured and organised output macro streams using __PRETTY_FUNCTION__
   INFO()    - info level messages, always displayed
   ERROR()   - error level messages, always displayed
   FATAL()   - error message that throws the process, always displayed 
   WARNING() - warning level messages, can be switched with the WARNINGLEVEL flag
   DEBUG()   - debug level messages, can be switched with the DEBUGLEVEL flag
   */

#include <iomanip>
#include <iostream>
#include <string>

#define WARNINGLEVEL 1
//#define DEBUGLEVEL 0
//#define TRACELEVEL 0
#define FCNNAMELENGTH 45

inline std::string trimmedString( std::string thing, const unsigned int& length = FCNNAMELENGTH )
{
  size_t pos2=0;
  do {
  pos2                            = thing.find( "AmpGen::" );
  if ( pos2 != std::string::npos ) thing = thing.replace( pos2, 8, "" );
  } while( pos2 != std::string::npos );

  pos2                                   = thing.find( "std::" );
  if ( pos2 != std::string::npos ) thing.replace( pos2, 5, "" );

  pos2                                   = thing.find( "virtual " );
  if ( pos2 != std::string::npos ) thing = thing.replace( pos2, 8, "" );

  size_t pos = thing.find( "(" );

  if ( pos != std::string::npos ) {
    return pos < length ? thing.substr( 0, pos ) : thing.substr( 0, length );
  }
  return thing.size() < length ? thing : thing.substr( 0, length ) + "...";
  if ( thing.size() < length ) return thing;
}

/// @ingroup msgService macro DEBUG
/// Used for printing verbose debugging messages, only if DEBUGLEVEL is defined.  
#ifdef DEBUGLEVEL
#define DEBUG( X )                                                                                                     \
  std::cout << "\033[2;32m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  DEBUG        "                                                                                       \
<< "\033[0m" << X << std::endl
#else
#define DEBUG( X )
#endif

/// @ingroup msgService macro INFO
/// Used for printing information messages, and will always be printed. 
#define INFO( X )                                                                                                      \
  std::cout << "\033[2;34m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  INFO         "                                                                                       \
<< "\033[0m" << X << std::endl

/// @ingroup msgService macro ERROR
/// Used for printing errors messages, and will always be printed. 
#define ERROR( X )                                                                                                     \
  std::cout << "\033[1;31m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  ERROR        "                                                                                                   \
<< "\033[0m" << X << std::endl

/// @ingroup msgService macro FATAL
/// Used for printing fatal errors messages, and will always be printed and will terminate the process afterwards.
#define FATAL( X )                                                                                                     \
  { std::cout << "\033[1;31m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  FATAL        "                                                                                                   \
<< "\033[0m" << X << std::endl;                                                                                         \
throw std::runtime_error( trimmedString( __PRETTY_FUNCTION__)+ " FATAL" ) ;}


/// @ingroup msgService macro FATAL
/// Used for printing warning messages, can be switched off using WARNINGLEVEL. These messages are often harmless, but sometimes not!
#ifdef WARNINGLEVEL
#define WARNING( X )                                                                                                   \
  std::cout << "\033[1;35m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  WARNING      "                                                                                       \
<< "\033[0m" << X << std::endl
#else
#define WARNING( X )
#endif

#ifdef TRACELEVEL
#define TRACE( X )                                                                                                     \
  std::cout << "\033[1;36m" << std::left << std::setw( FCNNAMELENGTH ) << trimmedString( __PRETTY_FUNCTION__ )         \
<< "  TRACE        "                                                                                       \
<< "\033[0m" << X << std::endl
#else
#define TRACE( X )
#endif

#endif
