#ifndef AMPGEN_SIMPLIFY_H
#define AMPGEN_SIMPLIFY_H 1
#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"

namespace AmpGen { 
  class NormalOrderedExpression { 
    public:
      struct Term {
        std::complex<double>                              m_prefactor; 
        std::vector<std::pair<Expression,std::string>>    m_terms;
        Expression                                        m_divisor; 
        std::string                                       m_expressionAsString; 
        bool                                              m_markForRemoval;
        void addExpression( const Expression& expression);
        Term( const Expression& expression ) ;
        operator Expression() ;
      };
      void addTerm( const Expression& expression );
      NormalOrderedExpression( const Expression& expression );
      void groupExpressions();

      operator Expression();
      std::vector<Expression> ExpandBrackets( const Expression& expression );
    private:
      std::vector<Term> m_terms; 
  };
  Expression Simplify(const Expression& expression );
}

#endif
