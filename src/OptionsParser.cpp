#include "AmpGen/OptionsParser.h"

#include <ctype.h>
#include <algorithm>
#include <utility>

#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MsgService.h"

using namespace AmpGen;

OptionsParser* OptionsParser::gOptionsParser = nullptr; 

OptionsParser::OptionsParser()
{
  m_keywords["Import"]             = [this](const auto& tokens){ if( tokens.size() != 2 ) return; this->import(expandGlobals(tokens[1] ) ) ; };
  m_keywords["ParticleProperties"] = [](const auto& tokens){ ParticlePropertiesList::getMutable()->addParticle( tokens ); };
  m_keywords["ParticlePropertiesList::Alias"] = [](const auto& tokens){ if( tokens.size() !=3 ) return; ParticlePropertiesList::getMutable()->makeAlias( tokens[1], tokens[2] ); };
}

OptionsParser* OptionsParser::getMe()
{
  if( gOptionsParser == nullptr ) gOptionsParser = new OptionsParser();
  return gOptionsParser;
}

void OptionsParser::setQuiet(){
  m_quiet = true; 
}

bool OptionsParser::ignoreThisLine( const std::string& line )
{
  if ( line.empty() ) return true;
  const char ignoreLinesStartingWith[] = {'*', '#', '\0'};
  for ( int i = 0; ignoreLinesStartingWith[i] != '\0'; i++ ) {
    if ( line[0] == ignoreLinesStartingWith[i] ) return true;
  }
  return false;
}

void OptionsParser::setCommandLineArgs( int argc, char** argv, const std::string& description )
{
  if( !m_quiet ) printSplash();
  int x = 0;
  while ( ++x < argc ) {
    std::string token = std::string( argv[x] );
    if ( token.find( "--" ) == std::string::npos ){
      import( token );   
      continue;
    }
    auto key        = token.substr( 2 );
    std::string val = "1";
    auto pos        = token.find( "=" );
    if ( pos != std::string::npos ) {
      key = token.substr( 2, pos - 2 );
      val = token.substr( pos + 1 );
    } else if ( x != argc - 1 ) {
      auto token2                                         = std::string( argv[x+1] );
      if ( token2.find( "--" ) == std::string::npos ){
        val = token2;
        x++;
      }
    }
    if( key == "help" ) m_printHelp = true; 
    addArg(key +" " + val); 
  }
  if( m_printHelp ){
    std::cout << bold_on << "Usage: " << bold_off << argv[0] << italic_on << " options_file1.opt options_file2.opt --key1=value1 --key2=value2 ..." << italic_off << std::endl; 
    if( description != "") std::cout << description << std::endl; 
    std::cout << bold_on << "Options: " << bold_off << std::endl; 
  }
}

void OptionsParser::import( const std::string& fName )
{
  if( !m_quiet) INFO( "Importing: " << fName );
  if ( !fileExists( fName ) ) {
    ERROR( "Cannot find file: " << fName );
    return;
  }
  int braceDepth = 0 ; 
  std::vector<std::string> currentTokens; 
  processFile( fName, [this, &currentTokens, &braceDepth]( auto& line ) {
    if ( this->ignoreThisLine( line ) ) return;
    auto tokens = this->makeParsedStrings( line, braceDepth );
    for ( auto& token : tokens ) currentTokens.push_back( token );
    if ( tokens.size() == 0 ) return;
    std::string key = currentTokens[0];
    if ( braceDepth != 0 ) return;
    if ( this->m_parsedLines.find( key ) != this->m_parsedLines.end() ) {
      WARNING( "Overwriting parameter: " << key );
    }
    auto isCurlyBrace = []( const std::string& o ) { return o == "{" || o == "}"; };
    currentTokens.erase( std::remove_if( currentTokens.begin(), currentTokens.end(), isCurlyBrace), currentTokens.end() );
    this->addArg( currentTokens );
    currentTokens.clear();
  } );
}

void OptionsParser::addArg( const std::string& arg )
{
  int bc = 0 ; 
  addArg( makeParsedStrings( arg, bc ) );
}

void OptionsParser::addArg( const std::vector<std::string>& tokens )
{
  auto& key = tokens[0];
  DEBUG("Adding arg with key: " << key );
  if( m_keywords.count(key) != 0 ) m_keywords[key]( tokens );
  else m_parsedLines[key] = tokens; 
}

std::vector<std::string> OptionsParser::makeParsedStrings( const std::string& line, int& braceDepth ) const
{
  std::string s = line;
  if ( s.empty() ) return {};
  s.push_back( ' ' );                         // makes sure we get last element
  s                                = " " + s; // makes things easier when we start with quotes.
  std::string::const_iterator prev = s.begin();
  bool prevBlank                   = true;
  bool insideQuotes                = false;
  bool ignore                      = false;
  std::vector<std::string> fillThisList;

  for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it ) {
    if ( !insideQuotes && *it == '#' ) break; /// indicates a comment, except in speech marks ///
    if ( ( ( !insideQuotes ) && std::isblank( *it ) ) || *it == '"' || ignore ) {
      if ( !prevBlank ) {
        std::string tmp_s( prev, it );
        fillThisList.push_back( tmp_s );
      }
      prevBlank = true;
    } else {
      if ( prevBlank ) prev = it;
      prevBlank             = false;
    }
    if ( *it == '"' ) insideQuotes = !insideQuotes;
    if ( *it == '{' && !insideQuotes ) braceDepth++;
    if ( *it == '}' && !insideQuotes ) braceDepth--;
  }
  if ( insideQuotes ) {
    WARNING( "Unbalanced quotes in string : " << line );
  }
  if ( braceDepth != 0 ) {
    DEBUG( "Opened braces in line: " << line );
  }
  return fillThisList;
}


bool                          OptionsParser::printHelp() { return getMe()->m_printHelp ; }   
void                          OptionsParser::setArgs( int argc, char** argv , const std::string& description){ getMe()->setCommandLineArgs(argc, argv, description); } 
void                          OptionsParser::setArg( const std::string& arg ){ getMe()->addArg( arg ); }
OptionsParser::iterator       OptionsParser::find( const std::string& name )  { return m_parsedLines.find( name ); }
OptionsParser::iterator       OptionsParser::begin()       { return m_parsedLines.begin(); }
OptionsParser::iterator       OptionsParser::end()         { return m_parsedLines.end(); }
OptionsParser::const_iterator OptionsParser::begin() const { return m_parsedLines.cbegin(); }
OptionsParser::const_iterator OptionsParser::end()   const { return m_parsedLines.cend(); }
