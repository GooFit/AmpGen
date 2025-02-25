#include "AmpGen/CompilerWrapper.h"

#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <ratio>
#include <dlfcn.h>
#include <unordered_map>
#include <map>
#include <utility>
#include <numeric>
#include <filesystem>

#include "AmpGen/NamedParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGenVersion.h"
#include "AmpGen/simd/utils.h"

using namespace AmpGen;
// #ifdef AMPGEN_CXX
// #pragma message "Using c++ compiler: " AMPGEN_CXX " for JIT"
// #pragma message "Using AMPGENROOT: "       AMPGENROOT
// #pragma message "Using AMPGENROOT_CMAKE: " AMPGENROOT_CMAKE 
// #else
// #pragma warning "No AMPGEN_CXX for JIT set"
// #endif 

CompilerWrapper::CompilerWrapper( const bool& verbose ) :
  m_verbose(verbose),
  m_cxx(getenv("AMPGEN_CXX") != nullptr ? std::string( getenv( "AMPGEN_CXX" ) ) : "")
{
  if ( m_cxx == "" ) {
#ifdef AMPGEN_CXX
    if( m_verbose ) 
      INFO( "Global variable AMPGEN_CXX is undefined, using compiler variable AMPGEN_CXX: " << AMPGEN_CXX );
    m_cxx = AMPGEN_CXX;
#else
    m_cxx = getenv("AMPGEN_CXX") != nullptr ? std::string( getenv( "AMPGEN_CXX" ) ) : "";
    if( m_cxx == "" ) ERROR( "No configured compiler; set global variable AMPGEN_CXX" );
#endif
  }
  else {
    if( m_verbose ) 
      INFO( "Global variable CXX is defined: " << m_cxx << "; AMPGEN_CXX is ignored" ); // to use AMPGEN_CXX, unset CXX
  }
}

void CompilerWrapper::generateSource( const CompiledExpressionBase& expression, const std::string& filename )
{
  std::ofstream output( filename );
  preamble(output);
  auto signature = expression.fcnSignature(); 
  if( signature.find("Complex")        != std::string::npos ){  output << "#include \"AmpGen/Complex.h\"\n" ; output << "using namespace AmpGen;\n";  }
  if( signature.find("AVX2d")        != std::string::npos )  output << "#include \"AmpGen/simd/avx2d_types.h\"\n; using namespace AmpGen::AVX2d;\n" ;
  else if( signature.find("AVX2f")   != std::string::npos )  output << "#include \"AmpGen/simd/avx2f_types.h\"\n; using namespace AmpGen::AVX2f;\n;" ;
  else if( signature.find("AVX512d") != std::string::npos )  output << "#include \"AmpGen/simd/avx512d_types.h\"\n; using namespace AmpGen::AVX512d;\n;" ;
  else if( signature.find("AVX512")  != std::string::npos )  output << "#include \"AmpGen/simd/avx512_types.h\"\n; using namespace AmpGen::AVX512;\n;" ;
  else if( signature.find("ARM128d") != std::string::npos )  output << "#include \"AmpGen/simd/arm128d_types.h\"\n; using namespace AmpGen::ARM128d;\n;" ;
  output << expression << std::endl; 
  output.close();
}

std::string CompilerWrapper::generateFilename()
{
  char buffer[] = "/tmp/libAmpGen-XXXXXX.cpp";
  int status    = mkstemps( buffer, 4);
  if ( status == -1 ) {
    ERROR( "Failed to generate temporary filename " << status );
  }
  std::remove( buffer );
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
  if ( name == "" ) name = generateFilename();  
  std::string cname = name +"_"+std::to_string(expression.hash())+".cpp";
  std::string oname = std::filesystem::path(cname).replace_extension(m_extension);
  if( NamedParameter<bool>("CompilerWrapper::ForceRebuild", false )  == false 
      && fileSize(oname) != -1 && expression.link( oname )) return true;
  if( print_all ) INFO("Generating source: " << cname );
  auto twall_begin  = std::chrono::high_resolution_clock::now();
  generateSource( expression, cname );
  compileSource(cname, oname);
  expression.link( oname );
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - twall_begin ).count();
  if( print_all ) INFO( expression.name() << " " << cname << " Compile time = " << tWall / 1000. << " [size = " << fileSize(cname)/1024 << ", " << fileSize(oname)/1024 << "] kB" );
  return true;
}


bool CompilerWrapper::compile( std::vector<CompiledExpressionBase*>& expressions, const std::string& fname, const std::map<std::string, std::string>& metadata_functions )
{
  bool print_all = m_verbose || NamedParameter<bool>("CompilerWrapper::Verbose",false);
  std::string cname = expandGlobals(fname);
  if ( cname == "" ) cname = generateFilename();
  std::string oname = std::filesystem::path(cname).replace_extension(m_extension);
  if( print_all ) INFO("Generating source: " << cname );
  auto twall_begin  = std::chrono::high_resolution_clock::now();
  std::ofstream output( cname );

  for ( const auto& include : m_includes ) output << "#include <" << include << ">\n";
  if( m_includePythonBindings ) output << "#include <cstring>\n"; 
  for ( const auto& expression : expressions ) output << *expression << "\n"; 
  for ( const auto& [k,v] : metadata_functions ) output << "extern \"C\" const char* " << k << "(){ return \"" << v << "\";}\n"; 

  output.close();
  compileSource( cname, oname );
  
  
  for( auto& expression : expressions ) expression->link( oname );
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - twall_begin ).count();
  if( print_all ) INFO( cname << " Compile time = " << tWall / 1000. << " [size = " << fileSize(cname)/1024 << ", " << fileSize(oname)/1024 << "] kB" );
  return true;
}

bool CompilerWrapper::isClang() const 
{
  return m_cxx.find("clang") != std::string::npos || m_cxx.find("llvm-g++") != std::string::npos;
}

std::string get_cpp_version()
{
  if( __cplusplus >= 201703L ) return "c++17";
  if( __cplusplus >= 201402L ) return "c++14";
  if( __cplusplus >= 201103L ) return "c++11";
  else return "";
}

void CompilerWrapper::compileSource( const std::string& fname, const std::string& oname )
{
  using namespace std::chrono_literals;
  std::vector<std::string> compile_flags = NamedParameter<std::string>("CompilerWrapper::Flags", {"-Ofast", "--std="+get_cpp_version()}); 
  compile_flags.push_back( std::string("-I") + AMPGENROOT) ;  
  #if INSTRUCTION_SET < 10
      compile_flags.push_back("-march=native");
    #if INSTRUCTION_SET == INSTRUCTION_SET_AVX2d 
      compile_flags.push_back("-mavx2");
      compile_flags.push_back("-DHAVE_AVX2_INSTRUCTIONS");
    #endif
  #endif
  if(std::string(AMPGEN_OPENMP_FLAGS) != ""){
    auto flags = split(AMPGEN_OPENMP_FLAGS, ' ');
    for( const auto& flag : flags ) compile_flags.push_back(flag); 
  }
  std::vector<const char*> argp = { m_cxx.c_str(), "-shared", "-rdynamic", "-fPIC"};  
  std::transform( compile_flags.begin(), compile_flags.end(), std::back_inserter(argp), [](const auto& flag ){return flag.c_str() ; } );
  
  if(isClang())
  {
    argp.push_back( "-Wno-return-type-c-linkage");
    #if __APPLE__
    argp.push_back("-lstdc++");
    argp.push_back("-isystem"); 
    #endif
    argp.push_back( "-march=native");
  }
  else {
    argp.push_back("-frename-registers");
  }
  auto tokens = split( std::string( AMPGEN_CXX_FLAGS ), ' ' );
  for( auto& token : tokens ) argp.push_back( token.c_str() );

  argp.push_back( fname.c_str() );
  argp.push_back( "-o");
  argp.push_back( oname.c_str() );

  if(NamedParameter<bool>("CompilerWrapper::Verbose", false)) {
    std::string result = std::accumulate(std::begin(argp), std::end(argp),
      std::string(),
      [](const std::string& a, const char* b){return a + " " + b;});
    INFO("Compiling: " << result);
  }
  argp.push_back( NULL );
  pid_t childPID = fork();
  if( childPID == 0 )
  {
    execv( argp[0], const_cast<char**>( &argp[0] ) );
    perror( "execv() in ../src/CompilerWrapper.cpp" ); // add "CompilerWrapper::Verbose true" to .opt file to diagnose
    exit( 0 );
  }
  int status = 0;   
  waitpid( childPID, &status, 0 );
}

void CompilerWrapper::preamble( std::ostream& os ) const 
{
  time_t now = time(0);
  char* dt = ctime(&now);
  const char* user = getenv("USER") == NULL ? "UNKNOWN" : getenv("USER");
  os << "/** Generated by " << user << " on " << dt ;
  os << " AmpGen v" << AMPGEN_MAJOR_VERSION << "." << AMPGEN_MINOR_VERSION << '\n'; 
  #if defined(__clang__)
    os << " clang v" << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
  #elif defined(__ICC) || defined(__INTEL_COMPILER)
    os << " icc " << __INTEL_COMPILER;
  #elif defined(__GNUC__) || defined(__GNUG__)
    os << " gcc v" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
  #endif
  os << '\n';
  os << "*/\n";

  for ( auto& include : m_includes ) os << "#include <" << include << ">\n";
}
