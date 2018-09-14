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
    std::string              m_cxx;
    std::string generateFilename();
  public:

    CompilerWrapper( const bool& verbose=false) ;

    void generateSource( const CompiledExpressionBase& expression, const std::string& fname);
    bool compile( CompiledExpressionBase& expression, const std::string& fname=""); 
    bool compile( std::vector<CompiledExpressionBase*>& expression, const std::string& fname=""); 
    void compileSource(const std::string& fname, const std::string& oname );

    void setVerbose() { m_verbose = true ; } 
  };
} // namespace AmpGen
#endif
