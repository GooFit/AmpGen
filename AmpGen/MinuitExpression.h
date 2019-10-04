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
    MinuitExpression( const std::vector<std::string>& tokens, MinuitParameterSet* mps );
    double getVal() const { return std::real( m_expression() ); }
    double mean()   const override { return getVal(); }
    operator double() const override { return getVal(); }
    ~MinuitExpression() override = default;
    bool isGood() const { return m_isGood; }
  
  private:
    Expression m_expression;
    bool m_isGood;
  };
} // namespace AmpGen

#endif
