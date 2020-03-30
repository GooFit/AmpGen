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
} // namespace AmpGen

#endif
