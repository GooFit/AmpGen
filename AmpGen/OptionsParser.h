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
  public:
    typedef std::map<std::string, std::vector<std::string>>::const_iterator const_iterator;
    typedef std::map<std::string, std::vector<std::string>>::iterator iterator;
    
    static OptionsParser* getMe();
    static bool printHelp();
    static void setArgs( int argc, char** argv, const std::string& description="" );
    static void setArg( const std::string& arg ); 
    void setQuiet(); 
    void addArg( const std::string& arg );
    void setCommandLineArgs( int argc, char** argv, const std::string& description =""); 
    void import( const std::string& fName );
    iterator find( const std::string& name );
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
  
  private:
    std::map<std::string,std::vector<std::string>>   m_parsedLines; 
    bool m_printHelp = {false};  
    bool m_quiet     = {false}; 
    static OptionsParser* gOptionsParser;
    
    OptionsParser() = default;
    bool ignoreThisLine( const std::string& line );
    void readStream( std::istream& is );
    std::vector<std::string> makeParsedStrings( const std::string& line, int& braceDepth ) const;
  };
} // namespace AmpGen
#endif
