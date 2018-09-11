#include "AmpGen/MinuitExpression.h"

#include "AmpGen/Expression.h"
#include "AmpGen/ExpressionParser.h"

namespace AmpGen
{
  class MinuitParameterSet;
} // namespace AmpGen

using namespace AmpGen;

MinuitExpression::MinuitExpression( const std::vector<std::string>& tokens, MinuitParameterSet* mps )
{
  setName( tokens[0] );
  std::string total_line = "";
  for ( size_t it = 2; it != tokens.size(); ++it ) total_line += tokens[it] + " ";
  ExpressionParser::getMe()->setMPS( mps );
  m_expression = ExpressionParser::Parse( total_line );
  m_isGood     = true;
}
