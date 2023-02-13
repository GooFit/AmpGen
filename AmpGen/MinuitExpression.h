#ifndef AMPGEN_MINUITEXPRESSION_H
#define AMPGEN_MINUITEXPRESSION_H

#include <complex>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/MinuitParameter.h"

namespace AmpGen
{
  class MinuitParameterSet;

  class MinuitExpression : public MinuitParameter
  {
  public:
    MinuitExpression(const std::vector<std::string>& tokens, MinuitParameterSet* mps );
    MinuitExpression(const std::string& name, const Expression& expression);
    double mean()   const override;
    complex_t getVal() const;
    Expression expression() const { return m_expression; }
    operator double() const override;
    ~MinuitExpression() override;  
  private:
    Expression m_expression;
  };
  
  class MinuitParameterLink : public IExpression {
    public:
      explicit MinuitParameterLink( MinuitParameter* param ) ;
      std::string to_string(const ASTResolver* resolver=nullptr) const override ;
      void resolve( ASTResolver& resolver ) const override ;
      complex_t operator()() const override ;
      operator Expression() const ;
      std::string name() const;
      const MinuitParameter& param() const;
    private: 
      MinuitParameter* m_parameter;
  };
} // namespace AmpGen

#endif
