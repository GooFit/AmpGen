#include "AmpGen/ASTResolver.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Spline.h"

using namespace AmpGen;

ASTResolver::ASTResolver(const std::map<std::string, size_t>& evtMap, 
    const MinuitParameterSet* mps ) : evtMap(evtMap), mps(mps), nParameters(0) {}

bool ASTResolver::hasSubExpressions() const { return subTrees.size() != 0; }
void ASTResolver::clearSubTrees() { subTrees.clear(); }

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
    parameter.m_address = res->second; 
    return;
  }
  else if( mps != nullptr ){
    auto it = mps->find(parameter.m_name);
    if( it != nullptr ){
      if( useCompileTimeConstants && it->iFixInit() == MinuitParameter::Flag::CompileTimeConstant ){
        parameter.m_defaultValue = it->mean();
        parameter.m_compileTimeConstant = true; 
      }
      else parameter.m_address = addCacheFunction<ParameterTransfer>( parameter.m_name, it );
      return;
    }
  }
  else if( useCompileTimeConstants ){
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
