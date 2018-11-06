#ifndef AMPGEN_OPTIONSPARSER_H
#define AMPGEN_OPTIONSPARSER_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <utility>

namespace AmpGen
{
  class OptionsParser
  {
  private:
    OptionsParser();
    static OptionsParser* gOptionsParser;
    bool ignoreThisLine( const std::string& line );
    void readStream( std::istream& is );
    void import( const std::string& fName );
    std::vector<std::string> makeParsedStrings( const std::string& line, int& braceDepth ) const;
  protected:
    std::map<std::string, std::vector< std::string> >   m_parsedLines; 
    bool m_printHelp;

  public:
    typedef std::map<std::string, std::vector<std::string>>::const_iterator const_iterator;
    typedef std::map<std::string, std::vector<std::string>>::iterator iterator;
    
    static OptionsParser* getMe();
    static bool printHelp() { return getMe()->m_printHelp ; }   
    static void setArgs( int argc, char** argv ){ getMe()->setCommandLineArgs(argc, argv ) ; } 
    static void setArg( const std::string& arg ){ getMe()->addArg( arg ); }
    void addArg( const std::string& arg );
    void setCommandLineArgs( int argc, char** argv );
     
    auto find( const std::string& name )  { return m_parsedLines.find( name ); }
    iterator begin() { return m_parsedLines.begin(); }
    iterator end() { return m_parsedLines.end(); }
    const_iterator begin() const { return m_parsedLines.cbegin(); }
    const_iterator end() const { return m_parsedLines.cend(); }
    
  };
} // namespace AmpGen
#endif
