#include "AmpGen/ASTResolver.h"

#include <algorithm>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Spline.h"
#include "AmpGen/Array.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/MinuitExpression.h"

using namespace AmpGen;

ASTResolver::ASTResolver(const std::map<std::string, unsigned> &evtMap, const MinuitParameterSet *mps) : m_evtMap(evtMap), m_mps(mps) {}

std::vector<std::pair<unsigned, Expression>> ASTResolver::getOrderedSubExpressions(const Expression &expression) {
  DEBUG("Calling solver on: " << expression);
  expression.resolve(*this);
  std::vector<std::pair<unsigned, Expression>> subexpressions;
  std::unordered_map<unsigned, size_t> used_functions;
  std::unordered_map<unsigned, Expression> sub_trees;
  do {
    sub_trees.clear();
    for(auto &[t, s] : m_tempTrees) {
      sub_trees[t] = cast<SubTree>(s).expression();
    }
    m_tempTrees.clear();
    for(auto &[key, sub_expression] : sub_trees) {
      sub_expression.resolve(*this);
      auto stack_pos = used_functions.find(key);
      if(stack_pos == used_functions.end()) {
        subexpressions.emplace_back(key, sub_expression);
        used_functions[key] = subexpressions.size() - 1;
        continue;
      }
      auto oldPos = stack_pos->second;
      auto it = subexpressions.begin() + oldPos;
      if(it == subexpressions.end() - 1) continue;
      std::rotate(it, it + 1, subexpressions.end());
      for(auto uf = used_functions.begin(); uf != used_functions.end(); ++uf) {
        if(uf->second >= oldPos) uf->second = uf->second - 1;
      }
      used_functions[key] = subexpressions.size() - 1;
    }
  } while(sub_trees.size() != 0);
  std::reverse(subexpressions.begin(), subexpressions.end());
  
  return subexpressions;
}

template <> void ASTResolver::resolve<SubTree>(const SubTree &subTree) {
  auto key = subTree.key();
  if(m_tempTrees.count(key) != 0)
    return; /// we already found this subtree, nothing to do
  else {
    subTree.expression().resolve(*this); // resolve the expression
    m_tempTrees[key] = subTree;         // and a pointer to this subtree for future consumption
  }
}

template <> void ASTResolver::resolve<Spline>(const Spline &spline) {
  if(m_resolvedParameters.count(&spline) != 0) return;
  auto address = addCacheFunction<SplineTransfer>(spline.m_name, spline.m_nKnots, spline.m_min, spline.m_max);
  addResolvedParameter(&spline, address);
  addResolvedParameter(spline.m_points.top().get(), address);
  auto splineTransfer = dynamic_cast<SplineTransfer *>(m_cacheFunctions[spline.m_name].get());
  if(m_mps == nullptr)
    ERROR("Fix me: Spline parameters must come from a ParameterSet");
  else {
    for(unsigned int i = 0; i < spline.m_nKnots; ++i) splineTransfer->set(i, m_mps->find(spline.m_name + "::" + std::to_string(i)));
  }
}

template <> void ASTResolver::resolve<Variable>(const Variable &parameter) {
  DEBUG("Resolving: " << parameter.name() << " " << parameter.isResolved());
  if(m_resolvedParameters.count(&parameter) != 0 || parameter.isResolved()) return;
  auto res = m_evtMap.find(parameter.name());
  if(res != m_evtMap.end()) {
    addResolvedParameter(&parameter, res->second, 1);
    return;
  } else if(m_mps != nullptr) {
    auto it = m_mps->find(parameter.name());
    if(it.parameter() != nullptr) {
      if(it->flag() == Flag::CompileTimeConstant) {
        addResolvedParameter(&parameter, "(" + std::to_string(it->mean()) + ")");
      } else
        addResolvedParameter(&parameter, addCacheFunction<ParameterTransfer>(parameter.name(), it));
      return;
    }
    DEBUG("Could not find parameter: " << parameter.name() << " amongst the MPS ");
  } else {
    addResolvedParameter(&parameter, std::to_string(parameter.defaultValue()));
    return;
  }
  auto address = addCacheFunction<CacheTransfer>(parameter.name(), parameter.defaultValue());
  addResolvedParameter(&parameter, address);
  if(!parameter.isResolved()) m_unresolvedVariables.push_back(parameter);
}

template <> void ASTResolver::resolve<ExpressionParameter>(const ExpressionParameter &parameter) {
  DEBUG("Resolving: " << parameter.name() << " " << m_resolvedParameters.count(&parameter));
  if(m_resolvedParameters.count(&parameter) != 0) return;
  if(m_mps == nullptr) return;
  auto it = m_mps->find(parameter.name());
  if(!it.isValid()) return;
  addResolvedParameter(&parameter, addCacheFunction<ParameterTransfer>(parameter.name(), it));
}

template <> void ASTResolver::resolve<LambdaExpression>(const LambdaExpression &obj) {
  if(m_resolvedParameters.count(&obj) != 0) return;
  addResolvedParameter(&obj, addCacheFunction<LambdaTransfer>("lxpr_" + std::to_string(m_nParameters), &obj));
}

std::map<std::string, std::shared_ptr<CacheTransfer>> ASTResolver::cacheFunctions() const { return m_cacheFunctions; }

void ASTResolver::addResolvedParameter(const IExpression *param, const std::string &thing) { m_resolvedParameters[param] = thing; }

void ASTResolver::addResolvedParameter(const IExpression *param, const size_t &address, const size_t &arg) {
  m_resolvedParameters[param] = "x" + std::to_string(arg) + "[" + std::to_string(address) + "]";
}

std::string ASTResolver::resolvedVariable(const IExpression *param) const {
  auto it = m_resolvedParameters.find(param);
  if(it != m_resolvedParameters.end())
    return it->second;
  else {
    ERROR("Parameter cannot be resolved: [" << param << "]: " << param->to_string());
    return "";
  }
}

void ASTResolver::clear() {
  m_resolvedParameters.clear();
  m_cacheFunctions.clear();
  m_parameterMapping.clear();
  m_tempTrees.clear();
  m_nParameters = 0;
}
