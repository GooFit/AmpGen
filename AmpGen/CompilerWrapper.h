#ifndef AMPGEN_COMPILERWRAPPER_H
#define AMPGEN_COMPILERWRAPPER_H

#include <dlfcn.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "AmpGen/CompiledExpressionBase.h"

namespace AmpGen
{
  class CompilerWrapper
  {
  private:
    std::vector<std::string> m_includes; 
    bool                     m_verbose;

    std::string generateFilename()
    {
      char buffer[] = "/tmp/libAmpGen-XXXXXX";
      int status    = mkstemp( buffer );
      if ( status == -1 ) {
        ERROR( "Failed to generate temporary filename " << status );
      }
      return buffer;
    }

  public:

    CompilerWrapper( const bool& verbose=false) : m_includes( {"array", "complex", "math.h", "vector"} ), m_verbose(verbose) {};

    void generateSource( const CompiledExpressionBase& expression, const std::string& fname);
    bool compile( CompiledExpressionBase& expression, const std::string& fname=""); 
    void setVerbose() { m_verbose = true ; } 
  };
} // namespace AmpGen
#endif
