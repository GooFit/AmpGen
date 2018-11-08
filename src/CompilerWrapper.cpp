#include "AmpGen/CompilerWrapper.h"

#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <bits/stdint-intn.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <ratio>
#include <dlfcn.h>
#include <unordered_map>
#include <map>
#include <utility>

#include "AmpGen/NamedParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CompiledExpressionBase.h"

using namespace AmpGen;

CompilerWrapper::CompilerWrapper( const bool& verbose ) : 
  m_includes( {"array", "complex", "math.h", "vector"} ),
  m_verbose(verbose)
{
  m_cxx = getenv("CXX") != nullptr ? std::string( getenv( "CXX" ) ) : "";
  if ( m_cxx == "" ) {
#ifdef AMPGEN_CXX
    if( m_verbose ) 
      INFO( "Using original compiler; set global variable CXX if another needed: " << AMPGEN_CXX );
    m_cxx = AMPGEN_CXX;
#else
    ERROR( "No configured compiler; set global variable CXX" );
#endif
  }
}

void CompilerWrapper::generateSource( const CompiledExpressionBase& expression, const std::string& filename )
{
  std::ofstream output( filename );
  for ( auto& include : m_includes ) output << "#include <" << include << ">\n";
  output << expression << std::endl; 
  output.close();
}

std::string CompilerWrapper::generateFilename()
{
  char buffer[] = "/tmp/libAmpGen-XXXXXX";
  int status    = mkstemp( buffer );
  if ( status == -1 ) {
    ERROR( "Failed to generate temporary filename " << status );
  }
  return buffer;
}

int64_t fileSize(const std::string& filename)
{
  struct stat stat_buf;
  int rc = stat(filename.c_str(), &stat_buf);
  return rc == 0 ? stat_buf.st_size : -1;
}


bool CompilerWrapper::compile( CompiledExpressionBase& expression, const std::string& fname )
{
  bool print_all = m_verbose || NamedParameter<bool>("CompilerWrapper::Verbose",false);

  std::string name = fname; 
  if ( name == "" ) 
    name = name == "" ? generateFilename() : expandGlobals( name );
  
  std::string cname = name +"_"+std::to_string(expression.hash())+".cpp";
  std::string oname = name +"_"+std::to_string(expression.hash())+".so";
  if( NamedParameter<bool>("CompilerWrapper::ForceRebuild", false )  == false 
      && fileSize(oname) != -1 && expression.link( oname )) return true;
  if( print_all ) INFO("Generating source: " << cname );
  auto twall_begin  = std::chrono::high_resolution_clock::now();
  generateSource( expression, cname );
  compileSource( cname, oname );
  expression.link( oname );
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - twall_begin ).count();
  if( print_all ) INFO( expression.name() << " " << cname << " Compile time = " << tWall / 1000. << " [size = " << fileSize(cname)/1024 << ", " << fileSize(oname)/1024 << "] kB" );
  return true;
}


bool CompilerWrapper::compile( std::vector<CompiledExpressionBase*>& expressions, const std::string& fname )
{
  bool print_all = m_verbose || NamedParameter<bool>("CompilerWrapper::Verbose",false);

  std::string name = fname;
  if ( name == "" ) name = name == "" ? generateFilename() : expandGlobals( name );
  std::string cname = name +".cpp";
  std::string oname = name +".so";
  if( print_all ) INFO("Generating source: " << cname );
  auto twall_begin  = std::chrono::high_resolution_clock::now();
  std::ofstream output( cname );
  for ( auto& include : m_includes ) output << "#include <" << include << ">\n";
  for ( auto& expression : expressions ) output << *expression << std::endl; 
  output.close();

  compileSource( cname, oname );
  for( auto& expression : expressions ) expression->link( oname );
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - twall_begin ).count();
  if( print_all ) INFO( cname << " Compile time = " << tWall / 1000. << " [size = " << fileSize(cname)/1024 << ", " << fileSize(oname)/1024 << "] kB" );
  return true;
}

void CompilerWrapper::compileSource( const std::string& fname, const std::string& oname )
{
  std::vector<pid_t> pids;
  pid_t childPID = 0; 
  using namespace std::chrono_literals;
  std::vector<const char*> argp = { m_cxx.c_str(), 
    "-Ofast", 
    "-shared", 
    "-rdynamic", 
    "--std=c++17",
    "-march=native",
    "-fPIC", 
    fname.c_str(), "-o", 
    oname.c_str(), NULL };
  childPID = vfork();
  if( childPID == 0 )
  {
    execv( argp[0], const_cast<char**>( &argp[0] ) );
    perror( "execl()" );
    exit( 0 );
  }
  int status = 0;   
  waitpid( childPID, &status, 0 );
}




