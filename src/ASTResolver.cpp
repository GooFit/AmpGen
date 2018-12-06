#include "AmpGen/ASTResolver.h"

#include <algorithm>
#include <cstdint>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Spline.h"
#include "AmpGen/Array.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

using namespace AmpGen;

ASTResolver::ASTResolver(const std::map<std::string, size_t>& evtMap, 
                         const MinuitParameterSet* mps ) : 
  evtMap(evtMap),
  mps(mps), 
  nParameters(0)
{
  enable_cuda = NamedParameter<bool>("enable_cuda",false);
  enable_compileTimeConstants = NamedParameter<bool>("enable_compileTimeConstants",false);
}

bool ASTResolver::hasSubExpressions() const 
{ 
  return subTrees.size() != 0; 
}

void ASTResolver::clearSubTrees() 
{ 
  subTrees.clear(); 
}

void ASTResolver::reduceSubTrees()
{
  subTrees.clear();
  for( auto& t : tempTrees ){
    auto expr = t.first->m_expression;
    uint64_t key = t.first->key(); 
    if( subTrees.count( key ) == 0 ){
      subTrees[ key ] = t.first->m_expression ;
    }
  }
  tempTrees.clear();
}

void ASTResolver::getOrderedSubExpressions( 
    Expression& expression, 
    std::vector< std::pair< uint64_t,Expression>>& dependentSubexpressions )
{
  expression.resolve( *this );
  std::map< uint64_t , size_t> used_functions;
  for( size_t i = 0 ; i < dependentSubexpressions.size(); ++i ) 
    used_functions[ dependentSubexpressions[i].first ] = i;
  do { 
    reduceSubTrees();
    for( auto& expression : subTrees ) {     
      expression.second.resolve( *this );
      auto stack_pos = used_functions.find( expression.first );
      if ( stack_pos == used_functions.end() ) {
        dependentSubexpressions.emplace_back( expression.first , expression.second );
        used_functions[expression.first] = dependentSubexpressions.size() - 1;
        continue;
      }
      unsigned int oldPos = stack_pos->second;
      auto it             = dependentSubexpressions.begin() + oldPos;
      if ( it == dependentSubexpressions.end() - 1 ) continue;

      std::rotate( it, it + 1, dependentSubexpressions.end() );

      for ( auto uf = used_functions.begin(); uf != used_functions.end(); ++uf ) {
        if ( uf->second >= oldPos ) uf->second = uf->second - 1;
      }
      used_functions[expression.first] = dependentSubexpressions.size() - 1;
    }
  } while ( hasSubExpressions() );
  std::reverse( dependentSubexpressions.begin(), dependentSubexpressions.end() );
}

template <> void ASTResolver::resolve<Parameter>( Parameter& parameter )
{
  auto res = evtMap.find(parameter.m_name);
  if( res != evtMap.end() ){
    if( enable_cuda ) {
      parameter.m_resolved = true;
      size_t t = res->second; 
      std::string it = ""; 
      if( t % 3 == 0 ) it = "x";
      if( t % 3 == 1 ) it = "y";
      if( t % 3 == 2 ) it = "z";
      int stg = t/3;
      if( stg == 0 ) parameter.m_name = "x1[i]."+it;
      else if( stg == 1 ) parameter.m_name = "x1[i+N]."+it;
      else parameter.m_name = "x1[i+"+std::to_string(stg)+"*N]." + it;
     // parameter.m_name = "x1[i+"+std::to_string(res->second)+"*N]";
    }
    else parameter.m_address = res->second; 
    return; 
  }
  else if( mps != nullptr ){
    auto it = mps->find(parameter.m_name);
    if( it != nullptr ){
      if( enable_compileTimeConstants && 
          it->iFixInit() == MinuitParameter::Flag::CompileTimeConstant ){
        parameter.m_defaultValue = it->mean();
        parameter.m_compileTimeConstant = true; 
      }
      else parameter.m_address = addCacheFunction<ParameterTransfer>( parameter.m_name, it );
      return;
    }
  }
  else if( enable_compileTimeConstants ){
    parameter.m_compileTimeConstant = true; 
    return; 
  }
  parameter.m_address = addCacheFunction<CacheTransfer>( parameter.m_name, parameter.m_defaultValue );
}

template <> void ASTResolver::resolve<SubTree>( SubTree& subTree )
{
  if( tempTrees.count( &subTree ) != 0 ) return;
  tempTrees[&subTree] = 1;
}

template <> void ASTResolver::resolve<Spline>( Spline& spline )
{
  spline.m_points.m_address = addCacheFunction<SplineTransfer>(spline.m_name,spline.m_nKnots,spline.m_min,spline.m_max);
  auto splineTransfer = dynamic_cast<SplineTransfer*>( m_cacheFunctions[spline.m_name].get() );
  if( mps == nullptr ) ERROR("Fix me!");
  if( spline.m_values.size() == 0 ){
    for( unsigned int i = 0 ; i < spline.m_nKnots; ++i ) 
      splineTransfer->set(i, mps->find(spline.m_name+"::"+std::to_string(i)) );
  }
  else for( unsigned int i = 0 ; i < spline.m_nKnots; ++i ) splineTransfer->set(i,spline.m_values[i]) ;
}
