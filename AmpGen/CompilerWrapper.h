#ifndef AMPGEN_COMPILERWRAPPER_H
#define AMPGEN_COMPILERWRAPPER_H

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace AmpGen
{
  class CompiledExpressionBase; 
  

  class CompilerWrapper
  {
  public:
    explicit CompilerWrapper( const bool& verbose=false);
    void generateSource( const CompiledExpressionBase& expression, const std::string& fname);
    bool compile( CompiledExpressionBase& expression, const std::string& fname=""); 
    bool compile( std::vector<CompiledExpressionBase*>& expression, const std::string& fname=""); 
    void compileSource(const std::string& fname, const std::string& oname );
    void setVerbose() { m_verbose = true ; } 
    void preamble(std::ostream& os ) const ; 
    void addHeader(const std::string& include ) { m_includes.push_back(include); } 
  
  private:
    std::vector<std::string> m_includes = {"complex", "cmath", "vector"}; 
    bool                     m_verbose;
    std::string              m_cxx;
    std::string generateFilename();
    bool isClang() const; 
    std::string              m_extension{"so"}; 
    
    class Cleaner {
      static Cleaner instance;
      Cleaner(Cleaner const&) = delete;             // Copy construct
      Cleaner(Cleaner&&) = delete;                  // Move construct
      Cleaner& operator=(Cleaner const&) = delete;  // Copy assign
      Cleaner& operator=(Cleaner &&) = delete;      // Move assign
      ~Cleaner(){
        std::cout << "Cleaning up objects..." << std::endl; 
      }
    };
  };
} // namespace AmpGen
#endif
