#ifndef AMPGEN_COMPILERWRAPPER_H
#define AMPGEN_COMPILERWRAPPER_H

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map> 
#include "AmpGen/Property.h" 

namespace AmpGen
{
  class CompiledExpressionBase; 
  std::string get_cpp_version(); 

  class CompilerWrapper
  {
    public:
      explicit CompilerWrapper();
      void generateSource( const CompiledExpressionBase& expression, const std::string& fname);
      bool compile( CompiledExpressionBase& expression, const std::string& fname=""); 
      bool compile( std::vector<CompiledExpressionBase*>& expression, 
          const std::string& fname="",
          const std::map<std::string, std::string>& metadata_functions = {} ); 
      void compileSource(const std::string& fname, const std::string& oname );
      void preamble(std::ostream& os ) const ; 
      void addHeader(const std::string& include ) { m_includes.push_back(include); } 
      void addPythonBindings(){ m_includePythonBindings = true ; }  
    private:
      using strings = std::vector<std::string>; 
      std::vector<std::string> m_includes = {"complex", "cmath", "vector"}; 
      bool                     m_includePythonBindings {true}; 
      std::string              m_cxx;
      std::string generateFilename();
      bool isClang() const; 
      std::string              m_extension{"so"};   
      Property<strings> m_compileFlags {this, "CompilerWrapper::Flags", {"-Ofast", "--std="+get_cpp_version()}}; 
      Property<bool>    m_verbose      {this, "CompilerWrapper::Verbose", false}; 
      Property<bool>    m_forceRebuild {this, "CompilerWrapper::ForceRebuild", false}; 
  };
} // namespace AmpGen
#endif
