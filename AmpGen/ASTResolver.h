#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <future>

#include "AmpGen/Expression.h"
#include "AmpGen/CacheTransfer.h"


namespace AmpGen {
  class MinuitParameter;
  class MinuitParameterSet;
  
  /** @class ASTResolver
   * Traverses trees in IExpression::resolveDependencies()
   * to keep track of the dependencies of the tree
   * such as sub-trees, external parameters and the event type map.
   */
  struct ASTResolver {
    std::map<std::string, std::shared_ptr<CacheTransfer>> cacheFunctions;
    std::map<std::string, size_t> evtMap;                /// event specification 
    std::map<std::string, std::string> parameterMapping; /// Mapping of parameters to compile parameters
    const MinuitParameterSet* mps;                       /// Set of MinuitParameters 
    std::map<SubTree*, int>       tempTrees;             /// temporary store of sub-trees for performing cse reduction 
    std::map<uint64_t, Expression> subTrees;             /// Unordered sub-trees 

    unsigned int nParameters; 
    bool useCompileTimeConstants;
    ASTResolver(const std::map<std::string, size_t>& evtMap = {} , 
        const MinuitParameterSet* mps = nullptr );
    bool hasSubExpressions() const;
    void clearSubTrees();

    unsigned int resolveParameter( const std::string& name, const double& defaultValue );
    void reduceSubTrees(); 
    void cleanup();
    template < class TYPE, class ...ARGS > size_t addCacheFunction( const std::string& name, ARGS&... args )
    {
      auto it = cacheFunctions.find(name);
      if( it != cacheFunctions.end() ) return it->second->address();
      auto cacheFunction = std::make_shared<TYPE>(nParameters, args... );
      cacheFunctions[name] = cacheFunction;
      nParameters += cacheFunction->size();
      return nParameters - cacheFunction->size();
    }
  };
}
