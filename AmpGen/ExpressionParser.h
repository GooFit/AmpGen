#ifndef AMPGEN_EXPRESSIONPARSER_H
#define AMPGEN_EXPRESSIONPARSER_H 1 

#include <complex>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/Types.h"
#include "AmpGen/enum.h"
#include "AmpGen/MsgService.h"

namespace AmpGen
{
  class ASTResolver;
  class MinuitParameter;
  class MinuitParameterSet;
  
  declare_enum(coordinateType, cartesian, polar)
  declare_enum(angType, deg, rad)
  
  template <typename arg_type> void set_tuple_from_expression( arg_type& arg, const Expression& expr )
  {
    if constexpr( std::is_same<arg_type, Expression>::value ){ arg = expr ; } 
    else if constexpr( std::is_same<arg_type, std::string>::value ){ arg = is<Parameter>(expr) ? cast<Parameter>(expr).name() : "ERROR"; }
    else if constexpr( std::is_same<arg_type, int>::value ){    arg = int(std::real(expr())); }
    else if constexpr( std::is_same<arg_type, double>::value ){ arg = double( std::real(expr())); } 
    else if constexpr( std::is_same<arg_type, std::complex<double>>::value ){ arg = expr(); }
    else FATAL( "Unrecognised type: " << type_string<arg_type> () ); 
  }

  template <typename... types> std::tuple< types... > unpack( const Expression& expression )
  {
    std::tuple<types...> rt; 
    unsigned int counter = 0; 
    if( ! is<ExpressionPack>( expression ) ) FATAL("Not a packed expression" );
    auto as_pack = cast<ExpressionPack>(expression).expressions(); 
    if( as_pack.size() != std::tuple_size<std::tuple<types...>>::value ) FATAL("Wrong number of expressions"); 
    for_each( rt, [&as_pack, &counter](auto& element)
    {
      set_tuple_from_expression(element, as_pack[counter] ); 
      counter++; 
    } );
    return rt; 
  }

  class ExpressionParser
  {
    public:
      typedef std::function<Expression( const Expression& )> unaryFCN;
      typedef std::function<Expression( const Expression&, const Expression& )> binaryFCN;
      
      void add_unary( const std::string& name, const unaryFCN& op ){ m_unaryFunctions[name] = op; }
      template <typename OP> void add_unary( const std::string& name )
      {
        add_unary(name, [](const auto& expression){ return OP(expression) ; } );
      }
      template<typename... args> void add_multi_arg( const std::string& name, const std::function<Expression(args...)>& fcn )
      {
        add_unary(name, [fcn](const auto& expression){ return std::apply(fcn, unpack<args...>(expression) ) ;} );
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
      bool                                           m_isCartesian = {true};
      double                                         m_sf          = {1};
    
  };
} // namespace AmpGen

#endif
