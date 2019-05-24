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
  m_isGood     = true;
  fix(); 
}
