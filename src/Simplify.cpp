#include "AmpGen/Simplify.h"

#include <stddef.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <functional>
#include <numeric>

#include "AmpGen/MsgService.h"

using namespace AmpGen;

void NormalOrderedExpression::Term::addExpression( const Expression& expression){
  if( is<IBinaryExpression>(expression) ){
    auto argAsBin = dynamic_cast<IBinaryExpression*>(expression.get());
    auto l = argAsBin->l(); 
    auto r = argAsBin->r(); 
    if( is<Product>(expression) ){
      addExpression(l);
      addExpression(r);
    }
    if( is<Sum>(expression) || is<Sub>(expression) ) 
      ERROR("I never should get here, brackets should be pre-expanded [" << expression << "] " );
    if( is<Divide>(expression) ){
      addExpression( l );
      m_divisor = m_divisor * r ;  
    }
  }
  else if( is<Constant>(expression) ) m_prefactor *= expression();
  else m_terms.emplace_back( expression, expression.to_string() );
}

NormalOrderedExpression::Term::Term( const Expression& expression ) : 
  m_prefactor(1), 
  m_divisor(1), 
  m_markForRemoval(false)
{
  addExpression(expression);
  std::sort( m_terms.begin(), 
             m_terms.end(), 
             [](auto& t1, auto& t2 ){ return t1.second > t2.second; } ) ;
  Expression t = std::accumulate( m_terms.begin(), m_terms.end(), Expression(1), [](auto& A, auto& B){ return A * B.first ; } ); 
 // std::multiplies<Expression>() );
  //Expression t = 1;
  //for( auto& f : m_terms ) t = t * f.first;
  m_expressionAsString = ( t / m_divisor).to_string();  
}

NormalOrderedExpression::Term::operator Expression() 
{ 
  Expression pf = m_prefactor;
  for( auto& t : m_terms ) pf = pf * t.first ;
  return pf / m_divisor; 
}

NormalOrderedExpression::NormalOrderedExpression(const Expression& expression, const bool& expandSubtrees ) : m_expandSubTrees(expandSubtrees)
{
  auto expanded = ExpandBrackets(expression);
  for( auto& t : expanded ){
    m_terms.emplace_back( t );
  }
  std::sort( m_terms.begin(), m_terms.end() , 
      [](auto& t1, auto&t2 ){ return t1.m_expressionAsString > t2.m_expressionAsString ; });
  groupExpressions();
}

void NormalOrderedExpression::groupExpressions(){
  for( size_t i= 0 ; i < m_terms.size() -1; ++i ){
    if( m_terms[i].m_markForRemoval ) continue; 
    for( size_t j = i + 1; j < m_terms.size(); ++j ){
      if( m_terms[i].m_expressionAsString != m_terms[j].m_expressionAsString ) break;
      m_terms[i].m_prefactor += m_terms[j].m_prefactor; 
      m_terms[j].m_markForRemoval = true; 
    }
  }
  m_terms.erase( std::remove_if( m_terms.begin(), m_terms.end(), [=]( auto& it ){ return it.m_markForRemoval ; }  ) , m_terms.end() );
}

NormalOrderedExpression::operator Expression(){
  Expression total = 0; 
  for( auto& t : m_terms ) total = total + t;  
  return total;
}

std::vector<Expression> NormalOrderedExpression::ExpandBrackets( const Expression& expression )
{
  if( is<IBinaryExpression>(expression) )
  {
    std::vector<Expression> rt; 
    if( is<Sum>(expression) ){
      auto left  = ExpandBrackets( cast<Sum>(expression).l() );
      auto right = ExpandBrackets( cast<Sum>(expression).r() );
      for( auto& l : left ) rt.emplace_back(l);
      for( auto& r : right ) rt.emplace_back(r);
      return rt; 
    }
    if( is<Sub>(expression) ){
      auto left  = ExpandBrackets( cast<Sub>(expression).l() );
      auto right = ExpandBrackets( cast<Sub>(expression).r() );
      for( auto& l : left ) rt.emplace_back(l);
      for( auto& r : right ) rt.emplace_back((-1)*r);
      return rt; 
    }
    if( is<Product>(expression) ){
      auto left  = ExpandBrackets( cast<Product>(expression).l() );
      auto right = ExpandBrackets( cast<Product>(expression).r() );
      for( auto& l : left ){
        for( auto& r : right ){
          rt.emplace_back(l*r);
        }
      }
      return rt; 
    }
    if( is<Divide>(expression) ){
      auto left  = ExpandBrackets( cast<Divide>(expression).l() );
      auto right = ExpandBrackets( cast<Divide>(expression).r() );
      Expression sum = 0;
      for( auto& r : right ) sum = sum + r; 
      for( auto& l : left ) rt.emplace_back(l/sum);
      return rt; 
    }
  }
  if( m_expandSubTrees && is<SubTree>(expression) ) 
    return ExpandBrackets( cast<SubTree>(expression).m_expression );
  if( is<Conj>(expression) ){
    auto terms = ExpandBrackets( cast<Conj>(expression).arg() );
    std::vector<Expression> rt;
    for(auto& it : terms) 
      rt.push_back( fcn::conj(it) );
    return rt;
  }
  return {expression}; 
}

Expression AmpGen::Simplify(const Expression& expression ){ return Expression( NormalOrderedExpression(expression) ) ; }
