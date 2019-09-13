#ifndef AMPGEN_COMPILEDEXPRESSIONBASE_H
#define AMPGEN_COMPILEDEXPRESSIONBASE_H
#include <memory.h>
#include <stddef.h>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <future>

#include "AmpGen/Expression.h"

namespace AmpGen
{
  class CacheTransfer;

  std::string programatic_name( std::string s );
  class ASTResolver;  
  class MinuitParameter;
  class MinuitParameterSet;

  /** @class CompiledExpressionBase
   *  Base class for compiled expressions, i.e. expressions that are (almost) ready to be evaluated.
   *  Handles (some) resolution and compilation behaviour, allows management of CompiledExpressions
   *  without explicitly referring to their return type, which is specified by template parameter
   *  in the implementation CompiledExpression
   */
  class CompiledExpressionBase
  {
  public:
    CompiledExpressionBase( const Expression& expression, 
                            const std::string& name,
                            const DebugSymbols& db=DebugSymbols(), 
                            const std::map<std::string,size_t>& evtMapping = {} );
    CompiledExpressionBase( const std::string& name );
    CompiledExpressionBase();

    void resolve(const MinuitParameterSet* mps = nullptr);
    void prepare();
    void compile(const std::string& fname=""); 
    void to_stream( std::ostream& stream ) const;
    unsigned int hash() const;
    std::string name() const;
    std::string progName() const;
    virtual bool link( void* handle              )          = 0;
    virtual bool link( const std::string& handle )          = 0;
    virtual void setExternal( const double& value, const unsigned int& address ) = 0;
    virtual void resizeExternalCache( const size_t& N ) = 0;
    virtual bool isReady() const               = 0;
    virtual std::string returnTypename() const = 0;
    virtual std::string fcnSignature()   const = 0;
    virtual std::string args()           const = 0;
    virtual void print() const                 = 0;
    virtual ~CompiledExpressionBase();
    virtual size_t returnTypeSize() const      = 0;    
  protected:
    Expression                                      m_obj;
    std::string                                     m_name;
    std::string                                     m_progName; 
    DebugSymbols                                    m_db;
    std::map<std::string, size_t>                   m_evtMap;
    std::vector<std::pair<uint64_t, Expression>>    m_dependentSubexpressions;
    std::vector<std::pair<uint64_t, Expression>>    m_debugSubexpressions; 
    std::vector<std::shared_ptr<CacheTransfer>>     m_cacheTransfers;
    std::shared_ptr<ASTResolver>                    m_resolver;
    
  private:
    void addDebug( std::ostream& stream ) const;
    void addDependentExpressions( std::ostream& stream, size_t& sizeOfStream ) const;
  };
  std::ostream& operator<<( std::ostream& os, const CompiledExpressionBase& expression );
}// namespace AmpGen
#endif
