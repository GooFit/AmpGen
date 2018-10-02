#ifndef AMPGEN_EXPRESSIONPARSER_H

#include <complex>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"


namespace AmpGen
{

  class MinuitParameterSet;
  struct ASTResolver;

  class MinuitParameterLink : public IExpression {
    public:
      MinuitParameterLink( MinuitParameter* param ) ;
      std::string to_string(const ASTResolver* resolver=nullptr) const override ;
      void resolve( ASTResolver& resolver ) override ;
      std::complex<double> operator()() const override ;
      operator Expression() const ;
      Expression clone() const { return MinuitParameterLink( m_parameter ) ; } 
    private: 
      MinuitParameter* m_parameter;
  };
  class ExpressionPack : public IExpression {
    public:
    ExpressionPack( const std::vector<Expression>& expressions ){ m_expressions = expressions ; }
    ExpressionPack( const Expression& A, const Expression& B );
    std::string to_string(const ASTResolver* resolver=nullptr) const override;
    void resolve( ASTResolver& resolver ) override;
    std::complex<double> operator()() const override;
    operator Expression() const ;
    Expression clone() const { 
      std::vector<Expression> cloned; 
      for( auto& expression : m_expressions ) cloned.push_back( expression.clone() );
      return ExpressionPack( cloned );
    }
    private:
    std::vector<Expression> m_expressions;
  };

  class ExpressionParser
  {
    private:
      typedef std::function<Expression( const Expression& )> unaryFCN;
      typedef std::function<Expression( const Expression&, const Expression& )> binaryFCN;

      MinuitParameterSet* m_mps;
      std::map<std::string, unaryFCN> m_unaryFunctions;
      std::vector<std::pair<std::string, binaryFCN>> m_binaryFunctions;

      void processBinaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions );

      void processUnaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions );

      Expression processEndPoint( const std::string& name );
      Expression parseTokens( const std::vector<std::string>& tokens );

      ExpressionParser();

    public:
      template <class OP>
        void add_unary( const std::string& name, OP op )
        {
          m_unaryFunctions[name] =
            std::function<Expression( const Expression& )>( [&op]( auto& expression ) { return op( expression ); } );
        }
      template <class OP>
        void add_unary( const std::string& name )
        {
          m_unaryFunctions[name] =
            std::function<Expression( const Expression& )>( []( auto& expression ) { return OP( expression ); } );
        }
      void add_binary( const std::string& name, binaryFCN op ) { m_binaryFunctions.emplace_back( name, op ); }

      void setMPS( MinuitParameterSet* mps ) { m_mps = mps; }
      static Expression Parse( const std::string& str );
      static ExpressionParser* gExpressionParser;
      static ExpressionParser* getMe()
      {
        if ( !gExpressionParser ) gExpressionParser = new ExpressionParser();
        return gExpressionParser;
      }
  };
} // namespace AmpGen

#endif
