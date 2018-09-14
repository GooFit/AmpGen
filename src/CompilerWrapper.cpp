#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/NamedParameter.h"

#include <dlfcn.h>
#include <stdio.h>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <sys/times.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

void CompilerWrapper::generateSource( const CompiledExpressionBase& expression, const std::string& filename )
{
  std::ofstream output( filename );
  for ( auto& include : m_includes ) output << "#include <" << include << ">\n";
  output << expression << std::endl; 
  output.close();
}

long GetFileSize(std::string filename)
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
  if( print_all ) INFO("Generating source: " << cname );
  auto twall_begin  = std::chrono::high_resolution_clock::now();
  generateSource( expression, cname );
  compileSource( cname, oname );
  expression.link( oname );
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - twall_begin ).count();
  if( print_all ) INFO( expression.name() << " " << cname << " Compile time = " << tWall / 1000. << " [size = " << GetFileSize(cname)/1024 << ", " << GetFileSize(oname)/1024 << "] kB" );
  return true;
}


bool CompilerWrapper::compile( std::vector<CompiledExpressionBase*>& expressions, const std::string& fname )
{
  bool print_all = true ; // m_verbose || NamedParameter<bool>("CompilerWrapper::Verbose",false);

  std::string name = fname;
  if ( name == "" ) 
    name = name == "" ? generateFilename() : expandGlobals( name );
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
  if( print_all ) INFO( cname << " Compile time = " << tWall / 1000. << " [size = " << GetFileSize(cname)/1024 << ", " << GetFileSize(oname)/1024 << "] kB" );
  return true;
}

void CompilerWrapper::compileSource( const std::string& fname, const std::string& oname )
{
  bool print_all = m_verbose || NamedParameter<bool>("CompilerWrapper::Verbose",false);
  const char* cxx            = getenv( "CXX" );
  if ( cxx == nullptr ) {
#ifdef AMPGEN_CXX
    if( print_all ) INFO( "Using original compiler; set global variable CXX if another needed: " << AMPGEN_CXX );
    cxx = AMPGEN_CXX;
#else
    ERROR( "No configured compiler; set global variable CXX" );
    return 0;
#endif
  }
  std::vector<pid_t> pids;
  pid_t childPID = 0; 
  using namespace std::chrono_literals;
  std::vector<const char*> argp = { cxx, 
    "-Ofast", 
    "-shared", 
    "-rdynamic", 
    "--std=c++14", 
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




