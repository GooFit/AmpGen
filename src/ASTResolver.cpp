#include "AmpGen/ASTResolver.h"
#include "AmpGen/MinuitParameterSet.h"

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
