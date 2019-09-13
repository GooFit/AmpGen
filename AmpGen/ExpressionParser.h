#ifndef AMPGEN_EXPRESSIONPARSER_H

#include <complex>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class ASTResolver;
  class MinuitParameter;
  class MinuitParameterSet;

  class MinuitParameterLink : public IExpression {
    public:
      MinuitParameterLink( MinuitParameter* param ) ;
      std::string to_string(const ASTResolver* resolver=nullptr) const override ;
      void resolve( ASTResolver& resolver ) const override ;
      complex_t operator()() const override ;
      operator Expression() const ;
      std::string name() const;
      const MinuitParameter& param() const;
    private: 
      MinuitParameter* m_parameter;
  };

  class ExpressionPack : public IExpression {
    public:
      ExpressionPack( const std::vector<Expression>& expressions ){ m_expressions = expressions ; }
      ExpressionPack( const Expression& A, const Expression& B );
      std::string to_string(const ASTResolver* resolver=nullptr) const override;
      void resolve( ASTResolver& resolver ) const override;
      complex_t operator()() const override;
      operator Expression() const ;
    private:
      std::vector<Expression> m_expressions;
  };

  class ExpressionParser
  {
    public:
      typedef std::function<Expression( const Expression& )> unaryFCN;
      typedef std::function<Expression( const Expression&, const Expression& )> binaryFCN;

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

      static Expression parse( const std::string& str, const MinuitParameterSet* mps=nullptr );
      static Expression parse(std::vector<std::string>::const_iterator begin,
                              std::vector<std::string>::const_iterator end  ,
                               const MinuitParameterSet* mps = nullptr ); 
      static ExpressionParser* gExpressionParser;
      static ExpressionParser* getMe()
      {
        if ( !gExpressionParser ) gExpressionParser = new ExpressionParser();
        return gExpressionParser;
      }
    private:
      ExpressionParser();

      Expression processEndPoint( const std::string& name, const MinuitParameterSet* mps = nullptr );
      Expression parseTokens( std::vector<std::string>::const_iterator begin, std::vector<std::string>::const_iterator end , const MinuitParameterSet* mps = nullptr );

      void processBinaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions );
      void processUnaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions );

      std::map<std::string, unaryFCN>                m_unaryFunctions;
      std::vector<std::pair<std::string, binaryFCN>> m_binaryFunctions;
  };
} // namespace AmpGen

#endif
