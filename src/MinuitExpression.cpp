#include "AmpGen/MinuitExpression.h"

#include <stddef.h>

#include "AmpGen/Expression.h"
#include "AmpGen/ExpressionParser.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

MinuitExpression::MinuitExpression( const std::vector<std::string>& tokens, MinuitParameterSet* mps )
{
  setName( tokens[0] );
  m_expression = ExpressionParser::parse(tokens.begin() + 2 , tokens.end() , mps );
  m_flag = Flag::Hide; 
}

MinuitExpression::MinuitExpression(const std::string& name, const Expression& expression)
{
  setName(name); 
  m_expression = expression; 
  m_flag = Flag::Hide; 
}

double MinuitExpression::mean() const { return std::real(getVal()); }
complex_t MinuitExpression::getVal() const { return m_expression(); } 
MinuitExpression::operator double() const { return std::real(getVal()); }
MinuitExpression::~MinuitExpression() = default;

