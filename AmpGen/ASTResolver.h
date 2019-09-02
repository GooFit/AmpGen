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
#include "AmpGen/CacheTransfer.h"
#include "AmpGen/Spline.h"

namespace AmpGen {
  class MinuitParameter;
  class MinuitParameterSet;
  class MinuitParameterLink;
  
  /** @ingroup ExpressionEngine class ASTResolver
      @brief (Internal) class to aide in the resolution of the dependencies of expression trees. 

      Traverses trees in IExpression::resolve() to keep track of the dependencies of the tree
      such as sub-trees, external parameters and the mapping between particle kinematics and the
      event vector.
  */ 
  class ASTResolver 
  { 
    public: 
      ASTResolver(const std::map<std::string, size_t>& evtMap = {} , const MinuitParameterSet* mps = nullptr );
      void cleanup();
      std::vector<std::pair<uint64_t,Expression>> getOrderedSubExpressions( const Expression& expression);

      template <class TYPE> void resolve( const TYPE& obj ){}
      template <class TYPE, class ...ARGS> size_t addCacheFunction( const std::string& name, const ARGS&... args )
      {
        auto it = m_cacheFunctions.find(name);
        if( it != m_cacheFunctions.end() ) return it->second->address();
        auto cacheFunction = std::make_shared<TYPE>(m_nParameters, args... );
        m_cacheFunctions[name] = cacheFunction;
        m_nParameters += cacheFunction->size();
        return m_nParameters - cacheFunction->size();
      }
      size_t nParams() const { return m_nParameters ; }       
      bool enableCuda() const { return m_enable_cuda ; }
      bool enableCompileConstants() const { return m_enable_compileTimeConstants ;} 
      std::map<std::string, std::shared_ptr<CacheTransfer>> cacheFunctions() const;
      void addResolvedParameter(const IExpression* param, const std::string& thing);
      void addResolvedParameter(const IExpression* param, const size_t& address, const size_t& arg=0);
      std::string resolvedParameter( const IExpression* param ) const; 

    private: 
      std::map<const IExpression*, std::string>             m_resolvedParameters;          /// Map of parameters that have been resolved
      std::map<std::string, std::shared_ptr<CacheTransfer>> m_cacheFunctions;              /// Container of functions for calculating function cache
      std::map<std::string, size_t>                         m_evtMap;                      /// Event specification 
      std::map<std::string, std::string>                    m_parameterMapping;            /// Mapping of parameters to compile parameters
      const MinuitParameterSet*                             m_mps;                         /// Set of MinuitParameters 
      std::map<const SubTree*, uint64_t>                    m_tempTrees;                   /// temporary store of sub-trees for performing cse reduction 
      unsigned int                                          m_nParameters;                 /// Number of parameters
      bool                                                  m_enable_cuda;                 /// flag to generate CUDA code <<experimental>>
      bool                                                  m_enable_compileTimeConstants; /// flag to enable compile time constants <<experimental>>
    
  };
  
  template <> void ASTResolver::resolve<Parameter>( const Parameter& obj );
  template <> void ASTResolver::resolve<SubTree>  ( const SubTree  & obj );
  template <> void ASTResolver::resolve<Spline>   ( const Spline   & obj );
  template <> void ASTResolver::resolve<MinuitParameterLink>( const MinuitParameterLink& obj );  
}
