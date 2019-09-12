#include "AmpGen/Utilities.h"
#include "AmpGenVersion.h"
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <utility>
#include <numeric>
#include "AmpGen/MsgService.h"

std::vector<std::string> AmpGen::vectorFromFile( const std::string& filename, const char ignoreLinesThatBeginWith )
{
  std::vector<std::string> output;
  std::string tmp;
  std::ifstream inFile( filename.c_str() );
  while ( inFile.good() ) {
    std::getline( inFile, tmp );
    if ( tmp.size() == 0 || tmp[0] == ignoreLinesThatBeginWith ) continue;
    output.push_back( tmp );
  }
  return output;
}

std::vector<std::string> AmpGen::split( const std::string& s, char delim, bool ignoreWhitespace )
{
  std::vector<std::string> elems;
  std::string item;
  std::stringstream ss( s );
  while ( std::getline( ss, item, delim ) ) {
    if ( !ignoreWhitespace || ( item != " " && item != "" && item != "\n" && item != "\t" ) ) elems.push_back( item );
  }
  return elems;
}

std::vector<std::string> AmpGen::split( const std::string& s, const std::vector<char>& delims )
{
  std::vector<std::string> elems;
  std::stringstream ss( s );
  std::string strDelim = std::accumulate( delims.begin(), delims.end(), std::string("") );
  std::string line;
  while ( std::getline( ss, line ) ) {
    std::size_t prev = 0, pos;
    while ( ( pos = line.find_first_of( strDelim, prev ) ) != std::string::npos ) {
      if ( pos > prev ) elems.push_back( line.substr( prev, pos - prev ) );
      prev = pos + 1;
    }
    if ( prev < line.length() ) elems.push_back( line.substr( prev, std::string::npos ) );
  }

  return elems;
}

void AmpGen::swapChars(std::string& arg, const char a, const char b)
{
  for( auto& c : arg )
  {
    if( c == a ) c = b;
    else if( c == b ) c = a;
  }
}

std::vector<size_t> AmpGen::findAll( const std::string& input, const std::string& ch )
{
  std::vector<size_t> output;
  size_t pos = 0;
  do {
    pos = input.find( ch, pos );
    if ( pos != std::string::npos ) {
      output.push_back( pos );
      pos++;
    }
  } while ( pos != std::string::npos );
  return output;
}

std::map<size_t, std::string> AmpGen::vecFindAll( const std::string& input, const std::vector<std::string>& vCh )
{
  std::map<size_t, std::string> output;
  for ( auto& ch : vCh ) {
    auto positions                          = findAll( input, ch );
    for ( auto& ip : positions ) output[ip] = ch;
  }
  return output;
}

std::string AmpGen::replaceAll( const std::string& input, const std::string& toReplace, const std::string& replaceWith )
{
  size_t pos         = 0;
  std::string output = input;
  if ( toReplace == replaceWith ) {
    ERROR( "This will lead to infinite loop and not do anything!" );
    return input;
  }
  size_t start_pos = 0;
  do {
    pos = output.find( toReplace, start_pos );
    if ( pos != std::string::npos ) {
      output.replace( pos, toReplace.length(), replaceWith );
      start_pos = pos + toReplace.length();
    }
  } while ( pos != std::string::npos );
  return output;
}

// extracts tree structures of the form X{Y,Z,A}
// where Y and Z and A are also tree elements, by finding
// the matching delimiter and the Z, A elements.

std::vector<std::string> AmpGen::getItems( const std::string& tree, 
                                           const std::vector<std::string>& brackets,
                                           const std::string& seperator )
{
  auto braces = vecFindAll( tree, brackets ); /// get a vector of positions of the brackets ///
  if ( braces.size() % 2 != 0 ) {
    ERROR( "Unmatched braces in expression: " << tree << " check string: " << braces.size() );
    return std::vector<std::string>();
  }

  if ( braces.size() == 0 ) return std::vector<std::string>( {tree} );
  std::vector<std::string> items = {tree.substr( 0, braces.begin()->first )};
  std::vector<std::pair<size_t, size_t>> matched_braces;
  for ( auto it = braces.begin(); it != braces.end(); ++it ) {
    const std::string& iType = it->second;
    if ( iType != brackets[0] ) continue;
    int diff = 1;
    for ( auto jt = it; jt != braces.end(); ++jt ) {
      if ( jt == it ) continue;
      const std::string& jType = jt->second; /// these are the sub-expressions ///
      diff += jType == brackets[1] ? -1 : +1;
      if ( diff != 0 ) continue;
      matched_braces.emplace_back( it->first, jt->first );
      break;
    }
  }
  auto commas = findAll( tree, seperator ); //// commas delimit separate decays;
  std::vector<std::string> daughterTrees;
  size_t begin_position = matched_braces.begin()->first + 1;
  for ( auto& comma : commas ) {
    auto itBrace = matched_braces.begin() + 1;
    for ( ; itBrace != matched_braces.end(); ++itBrace ) {
      if ( comma > itBrace->first && comma < itBrace->second ) break;
    }
    if ( itBrace == matched_braces.end() ) {
      items.push_back( tree.substr( begin_position, comma - begin_position ) );
      begin_position = comma + 1;
    }
  }
  items.push_back( tree.substr( begin_position, matched_braces.begin()->second - begin_position ) );
  return items;
}

size_t AmpGen::find_next_of( const std::string& input, const std::vector<std::string>& patterns, const size_t& begin )
{
  size_t minPos = std::string::npos;
  for ( auto& pattern : patterns ) {
    size_t pos                 = input.find( pattern, begin );
    if ( pos < minPos ) minPos = pos;
  }
  return minPos;
}

unsigned int AmpGen::editDistance( const std::string& s1, const std::string& s2 )
{
  const std::size_t len1 = s1.size(), len2 = s2.size();
  std::vector<std::vector<unsigned int>> d( len1 + 1, std::vector<unsigned int>( len2 + 1 ) );

  d[0][0] = 0;
  for ( unsigned int i = 1; i <= len1; ++i ) d[i][0] = i;
  for ( unsigned int i = 1; i <= len2; ++i ) d[0][i] = i;

  for ( unsigned int i = 1; i <= len1; ++i )
    for ( unsigned int j = 1; j <= len2; ++j )
      d[i][j] = std::min( {d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + ( s1[i - 1] == s2[j - 1] ? 0 : 1 )} );
  return d[len1][len2];
}

std::string AmpGen::round( const double& number, const unsigned int& nsf )
{
  double value = std::round(number * pow(10, nsf)) / pow(10, nsf);
  auto r = mysprintf( ( "%." + std::to_string(nsf) + "f" ).c_str(), value );
  return r.substr(0,r.size()-1);
}

std::string AmpGen::numberWithError( const double& number, const double& error, const unsigned int& nDigits )
{
  return round( number, nDigits ) + "\\pm" + round( error * pow( 10, nDigits ), 0 );
}

bool AmpGen::stringMatchesWildcard( const std::string& input, const std::string& wildcard_string,
                            const char wildcard_character )
{
  auto pos = wildcard_string.find( wildcard_character ); /// TEST_foobar -> *_foobar
  if ( wildcard_string.size() == 1 && wildcard_string[0] == wildcard_character ) {
    DEBUG( "Returning true" );
    return true;
  }
  if ( pos == std::string::npos ) {
    DEBUG( "Returning " << input << " = " << wildcard_string << " ?" );
    return input == wildcard_string;
  }
  if ( pos == wildcard_string.size() - 1 ) {
    DEBUG( "Returning " << input << " contains " << wildcard_string );
    return input.find( wildcard_string.substr( 0, wildcard_string.size() - 1 ) ) == 0;
  } else {
    const std::string pattern1 = wildcard_string.substr( 0, pos + 1 );
    const std::string pattern2 = wildcard_string.substr( pos + 1 );
    DEBUG( "Matching " << pattern1 << " to " << input );
    bool match1 = stringMatchesWildcard( input, pattern1, wildcard_character );
    if ( !match1 ) return false;
    auto pos2             = pattern2.find( wildcard_character );
    auto posInInputString = input.find( pattern2.substr( 0, pos2 ) );
    if ( posInInputString == std::string::npos ) return false;
    return stringMatchesWildcard( input.substr( posInInputString ), pattern2, wildcard_character );
  }
  return false;
}

unsigned int AmpGen::FNV1a_hash( const std::string& toHash )
{ 
  unsigned int hash = 2166136261;
  for ( auto& c : toHash ) {
    hash ^= c;
    hash *= 16777619;
  }
  return hash;
}

std::ostream& AmpGen::bold_on(std::ostream& os){ return os << "\033[1m"; }

std::ostream& AmpGen::bold_off(std::ostream& os){ return os << "\033[0m";}

std::ostream& AmpGen::italic_on(std::ostream& os){ return os << "\033[3m"; }

std::ostream& AmpGen::italic_off(std::ostream& os){ return os << "\033[0m";}

void AmpGen::printReleaseNotes( const std::string& fname )
{
  bool printLines = false; 
  auto lines = vectorFromFile( fname ); 
  for( auto& line : lines )
  {
    if( line.find("!=") == 0 && ! printLines ){ printLines = true; 
      std::cout << bold_on << " Release notes (" << split( line, ' ' )[2] << ")" << bold_off << std::endl; 
      continue; 
    }
    if( ! printLines ) continue; 
    if( line.find(" -") == 0 ) std::cout << line << std::endl; 
    if( line.find("!=") == 0 ) break;  
  }
  std::cout << std::endl; 
}


void AmpGen::printSplash()
{
  std::cout << "\n\033[2;31m";
  std::cout << "    █████╗ ███╗   ███╗██████╗  ██████╗ ███████╗███╗   ██╗" << std::endl;
  std::cout << "   ██╔══██╗████╗ ████║██╔══██╗██╔════╝ ██╔════╝████╗  ██║" << std::endl;
  std::cout << "   ███████║██╔████╔██║██████╔╝██║  ███╗█████╗  ██╔██╗ ██║" << std::endl;
  std::cout << "   ██╔══██║██║╚██╔╝██║██╔═══╝ ██║   ██║██╔══╝  ██║╚██╗██║" << std::endl;
  std::cout << "   ██║  ██║██║ ╚═╝ ██║██║     ╚██████╔╝███████╗██║ ╚████║" << std::endl;
  std::cout << "   ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝      ╚═════╝ ╚══════╝╚═╝  ╚═══╝" << std::endl;
  std::cout << "\033[0m\n";
  std::cout << bold_on << 
               "                      version " << AMPGEN_MAJOR_VERSION << "." << AMPGEN_MINOR_VERSION << std::endl; 
  std::cout << "        build: " ;
  #if defined(__clang__)
    std::cout << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
  #elif defined(__llvm__)
    std::cout<< "llvm ";
  #elif defined(__ICC) || defined(__INTEL_COMPILER)
      std::cout << "icc " << __INTEL_COMPILER;
  #elif defined(__GNUC__) || defined(__GNUG__)
    std::cout << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
  #endif
  std::cout << "  " << __DATE__ << " " << __TIME__ << bold_off << "\n\n";

  char* AmpGenRoot = getenv("AMPGENROOT");
  if( AmpGenRoot != nullptr ) printReleaseNotes( std::string(AmpGenRoot) + "/doc/release.notes"); 
}

bool AmpGen::fileExists( const std::string& name )
{
  struct stat buffer;
  return ( stat( name.c_str(), &buffer ) == 0 );
}

std::string AmpGen::rtrim( std::string s )
{
  s.erase( std::find_if( s.rbegin(), s.rend(), [](unsigned char c){ return !std::isspace(c); }).base(), s.end() );
  return s;
}

// trim from both ends
std::string AmpGen::trim( std::string s ) { return ltrim( rtrim( s ) ); }

std::string AmpGen::ltrim( std::string s )
{
  s.erase( s.begin(), std::find_if( s.begin(), s.end(), [](unsigned char c){ return !std::isspace(c); }));
  return s;
}

std::string AmpGen::expandGlobals( std::string path )
{
  size_t pos;
  do {
    pos = path.find( "$" );
    if ( pos == std::string::npos ) break;
    size_t end_pos; // = std::string::npos;
    std::string variable_name;
    if ( path[pos + 1] == '{' ) {
      end_pos       = path.find( "}", pos );
      variable_name = path.substr( pos + 2, end_pos - pos - 2 );
      end_pos       = end_pos + 1;
    } else {
      end_pos       = find_next_of( path, {".", "/"}, pos );
      variable_name = path.substr( pos + 1, end_pos - pos - 1 );
    }
    const char* global_var = getenv( variable_name.c_str() );
    if ( global_var == nullptr ) {
      ERROR( "variable " << variable_name << " not found" );
      break;
    }
    std::string old_path = path;
    size_t len           = end_pos == std::string::npos ? path.length() - pos + 1 : end_pos - pos + 1;
    path                 = path.replace( pos, len, global_var );
    DEBUG( old_path << " -> " << path );
  } while ( pos != std::string::npos );

  return path;
}

bool AmpGen::isDir( const std::string& pathname )
{
  struct stat sb;
  return stat( pathname.c_str(), &sb ) == 0 && S_ISDIR( sb.st_mode );
}

std::vector<std::string> AmpGen::getListOfFiles( const std::string& directory, const std::string& patternString )
{
  std::string expanded_path = expandGlobals( directory );
  std::vector<std::string> files;
  std::vector<std::string> top_paths;

  auto wildcard_pos = expanded_path.find( "*" );
  if ( wildcard_pos != std::string::npos ) {
    auto end                   = expanded_path.find( "/", wildcard_pos );
    auto begin                 = expanded_path.find_last_of( "/", wildcard_pos );
    std::string sub_path       = expanded_path.substr( 0, begin );
    std::string matchingString = expanded_path.substr( begin + 1, end - begin - 1 );
    auto theseFiles            = getListOfFiles( sub_path, matchingString );
    for ( auto& file : theseFiles ) {
      if ( isDir( file ) )
        top_paths.push_back( file );
      else
        files.push_back( file );
    }
  } else
    top_paths.push_back( expanded_path );
  struct dirent* ent;
  for ( auto& top_path : top_paths ) {
    DIR* dir;
    if ( ( dir = opendir( top_path.c_str() ) ) != nullptr ) {
      /* print all the files and directories within directory */
      while ( ( ent = readdir( dir ) ) != nullptr ) {
        // printf ("%s\n", ent->d_name);
        std::string name = ent->d_name;
        if ( name == ".." || name == "." ) continue;
        if ( patternString == "" || stringMatchesWildcard( name, patternString ) )
          files.push_back( top_path + "/" + name );
      }
      closedir( dir );
    } else {
      /* could not open directory */
      perror( "" );
      return files;
    }
  }
  return files;
}

template <>
double AmpGen::lexical_cast<double>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtod( word.c_str(), &p );
  status &= *p == 0;
  return number;
}

template <>
int AmpGen::lexical_cast<int>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtol( word.c_str(), &p, 10 );
  status &= *p == 0;
  return number;
}
template <>
unsigned int AmpGen::lexical_cast<unsigned int>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtoul( word.c_str(), &p, 10 );
  status &= *p == 0;
  return number;
}

template <>
std::string AmpGen::lexical_cast<std::string>( const std::string& word, bool& status )
{
  return word;
}

template <>
float AmpGen::lexical_cast<float>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtof( word.c_str(), &p );
  status &= *p == 0;
  return number;
}

template <>
bool AmpGen::lexical_cast<bool>( const std::string& word, bool& status )
{
  bool value = false;
  if ( word == "1" || word == "true" ) {
    value = true;
  } else if ( word == "0" || word == "false" ) {
    value = false;
  } else
    status = false;

  return value;
}

template <>
long int AmpGen::lexical_cast<long int>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtol( word.c_str(), &p, 10 );
  status &= *p == 0;
  return number;
}

template <>
unsigned long int AmpGen::lexical_cast<unsigned long int>( const std::string& word, bool& status )
{
  char* p;
  auto number = strtoul( word.c_str(), &p, 10 );
  status &= *p == 0;
  return number;
}

std::pair<size_t,int> AmpGen::minSwaps(const std::vector<size_t>& indices, const std::vector<int>& exchangeParities) 
{ 
  std::vector<std::pair<size_t, size_t>> arrPos;
  for(size_t i = 0; i < indices.size(); i++) arrPos.emplace_back( indices[i], i ); 
  std::sort(arrPos.begin(), arrPos.end() ); 
  std::vector<bool> visited(indices.size(), false); 
  size_t ans = 0; 
  int netExchangeParity = +1;
  for (size_t i = 0; i < indices.size(); i++) 
  { 
    if (visited[i] || arrPos[i].second == i)  continue; 
    size_t cycle_size = 0; 
    size_t j = i; 
    while (!visited[j]) 
    { 
      visited[j] = true; 
      j = arrPos[j].second; 
      cycle_size++; 
    } 
    if( cycle_size > 0 ){
      netExchangeParity *= exchangeParities[i]; 
      ans += (cycle_size - 1); 
    }
  } 
  return {ans,netExchangeParity}; 
}
