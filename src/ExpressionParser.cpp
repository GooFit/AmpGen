#include "AmpGen/ExpressionParser.h"

#include <stddef.h>
#include <cmath>
#include <ostream>
#include <functional>

#include <TRandom3.h>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitExpression.h"
#include "AmpGen/ASTResolver.h"
#include "AmpGen/enum.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Tensor.h"

using namespace AmpGen;
using namespace std::complex_literals; 

DEFINE_CAST( MinuitParameterLink )

namespace AmpGen 
{
  complete_enum(coordinateType, cartesian, polar)
  complete_enum(angType, deg, rad)
}

void ExpressionParser::processBinaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions )
{
  for ( auto& fcn : m_binaryFunctions ) {
    for ( int pos = 0; pos < int(opCodes.size()); ++pos ) {
      if ( opCodes[pos] != fcn.first ) continue;
      expressions[pos] = fcn.second( expressions[pos], expressions[pos + 1] );
      expressions.erase( expressions.begin() + pos + 1 );
      opCodes.erase( opCodes.begin() + pos );
      pos--;
    }
  }
}

void ExpressionParser::processUnaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions )
{
  for ( int pos = 0; pos < int(opCodes.size()); ++pos ) {
    auto fcn = m_unaryFunctions.find( opCodes[pos] );
    DEBUG( " op = " << opCodes[pos] << " pos = " << pos );
    if ( fcn == m_unaryFunctions.end() ) continue;
    expressions[pos] = fcn->second( expressions[pos] );
    opCodes.erase( opCodes.begin() + pos );
    pos--;
  }
}

using iterator = std::vector<std::string>::const_iterator; 

template <char open= '(', char close = ')' > iterator findMatchingBracket(iterator begin, iterator end )
{
  int openedBrackets = 1;
  if( begin + 1 >= end ) return end; 
  for( auto it = begin+1; it != end; ++it )
  {
    openedBrackets += ( it->at(0) == open) -( it->at(0) == close );
    if( openedBrackets == 0 ) return it;
  }
  for( auto it = begin; it != end; ++it )
  {
    std::cout << *it << " ";
  }
  std::cout << std::endl; 
  FATAL("Unmatched brace in expression, " << open << " " << close);
  return end; 
}

/// divides strings by a delimiter, allowing for bracketing rules 
std::vector< std::pair<iterator, iterator> > split_it( iterator begin, iterator end, const std::string& delimiter)
{
  std::vector<std::pair<iterator, iterator>> rt; 
  rt.emplace_back(begin, begin);
  int b1l = 0;
  int b2l = 0;
  for( auto it = begin; it != end-1; ++it )
  {
    b1l += ( it->at(0) == '(') - (it->at(0) == ')');
    b2l += ( it->at(0) == '{') - (it->at(0) == '}');
    if( *it == delimiter && b1l == 0 && b2l == 0 ) 
    {
      rt.rbegin()->second = it; 
      rt.emplace_back( it+1, it+1);
    }
  }
  rt.rbegin()->second = end; 
  return rt; 
} 

std::string merge( iterator begin, iterator end)
{
  std::string total = "";
  for( auto it = begin; it != end; ++it ) total += *it + " ";
  return total;
}

Expression ExpressionParser::parseTokens(std::vector<std::string>::const_iterator begin,
                                         std::vector<std::string>::const_iterator end,
                                         const MinuitParameterSet* mps )
{
  std::vector<std::string> opCodes;
  std::vector<Expression> expressions; 
  for( auto it = begin; it != end; ++it )
  {
    if( it->at(0) == '(' ){
      auto begin_2 = it;
      auto end_2   = findMatchingBracket<'(',')'>(it, end);
      expressions.emplace_back( parseTokens(begin_2+1, end_2, mps) );
      it = end_2;
    }
    else if ( it->at(0) == '{' )
    {
      auto begin_2 = it; 
      auto end_2   = findMatchingBracket<'{','}'>(it, end);
      auto divided = split_it( begin_2+1, end_2, ",");
      Tensor v( std::vector<unsigned>{static_cast<unsigned>(divided.size())} );
      for( unsigned int i = 0 ; i != divided.size(); ++i )
      {
        v[i] = parseTokens( divided[i].first, divided[i].second, mps); 
      }
      DEBUG( "Adding tensor: " << v.to_string() );
      expressions.emplace_back(TensorExpression(v));
      it = end_2;
    }
    else {
      auto f = std::find_if(m_binaryFunctions.begin(), m_binaryFunctions.end(), [it](auto& jt){ return jt.first == *it; } );
      if( f != m_binaryFunctions.end() || m_unaryFunctions.count(*it) )
        opCodes.push_back(*it);
      else expressions.push_back( processEndPoint( *it , mps ) );
    }
  }
  DEBUG( "nExpressions = " << expressions.size() << " nOpCodes = " << opCodes.size() << " " << vectorToString( opCodes, " " ) );
  processUnaryOperators( opCodes, expressions );
  processBinaryOperators( opCodes, expressions );
  if( expressions.size() != 1 ){
    ERROR("Could not process expression: n = " << expressions.size() << " " << merge(begin, end) );
  }
  return expressions[0];
}

ExpressionParser::ExpressionParser()
{
  add_unary<Sin>( "sin" ); /// "function" operator ordering is irrelevant
  add_unary<Cos>( "cos" );
  add_unary<Tan>( "tan" );

  add_unary<ACos>( "acos" );
  add_unary<ATan>( "atan" );
  add_unary<ASin>( "asin" );

  add_unary<Sqrt>( "sqrt" );
  add_unary<Exp>( "exp" );
  add_unary<Log>( "log" );
  add_unary<Real>("real");
  add_unary<Imag>("imag");
  add_unary<Abs>("abs");
  
  add_multi_arg("atan2"         , std::function<Expression(Expression, Expression)>( [](Expression A, Expression B){ return fcn::atan2(A, B); }) );  
  add_multi_arg("Ternary"       , std::function<Expression(Expression, Expression, Expression)>( [](Expression A, Expression B, Expression C){ return Ternary(A, B, C) ; } ) ); 
  add_multi_arg("Secret", std::function<Expression(std::string, double, double)>( [](std::string key, double min, double max){ return TRandom3( FNV1a_hash(key) ).Uniform(min, max); } ) );

  /// the order of operations matters here
  add_binary( "^" , [](const auto& A, const auto& B ) { return fcn::pow(A,B); } );
  add_binary( "/" , [](const auto& A, const auto& B ) { return A / B; } ); 
  add_binary( "*" , [](const auto& A, const auto& B ) { return A * B; } );
  add_binary( "-" , [](const auto& A, const auto& B ) { return A - B; } );
  add_binary( "+" , [](const auto& A, const auto& B ) { return A + B; } );
  add_binary( ">" , [](const auto& A, const auto& B ) { return A > B; } );
  add_binary( "<" , [](const auto& A, const auto& B ) { return A < B; } );
  add_binary( "&&", [](const auto& A, const auto& B ) { return A && B; } );
  add_binary( "," , [](const auto& A, const auto& B ) { return ExpressionPack( A, B ); } );
 
  coordinateType coord = NamedParameter<coordinateType>("CouplingConstant::Coordinates", coordinateType::cartesian);
  angType degOrRad     = NamedParameter<angType>("CouplingConstant::AngularUnits", angType::rad);
  m_isCartesian = true; 
  if( coord == coordinateType::polar ) m_isCartesian = false; 
  else if ( coord != coordinateType::cartesian){
    FATAL("Coordinates for coupling constants must be either cartesian or polar");
  } 
  if ( degOrRad == angType::deg) m_sf = M_PI / 180; 
  else if ( degOrRad != angType::rad){
    FATAL("CouplingConstant::AngularUnits must be either rad or deg");
  } 
}

Expression ExpressionParser::parse( 
    std::vector<std::string>::const_iterator begin,
    std::vector<std::string>::const_iterator end,
    const MinuitParameterSet* mps ) { 
  auto parsed = getMe()->parseTokens(begin, end, mps ); 
  ASTResolver resolver( {}, mps ); 
  parsed.resolve( resolver ); 
  if( resolver.unresolvedParameters().size() != 0 )
  {
    auto expr = vectorToString( begin, end, " ", [](const auto& str){ return str; } ); 
    ERROR("Expression: " << expr << " contains unknown parameters:");
    for( const auto& parameter : resolver.unresolvedParameters() ) ERROR ( parameter.name() );
  }
  return parsed; 
}

Expression ExpressionParser::parse( 
    const std::string& expr,
    const MinuitParameterSet* mps ) { 
  auto tokens = split( expr , ' ' );
  return getMe()->parseTokens(tokens.cbegin(), tokens.cend() , mps ); 
}

ExpressionParser* ExpressionParser::gExpressionParser = nullptr;

Expression ExpressionParser::processEndPoint( const std::string& name, const MinuitParameterSet* mps )
{
  bool status  = true;
  double value = lexical_cast<double>( name, status );
  if ( status == true ) return value;
  if ( name == "PI" || name == "pi" || name == "M_PI" ) return M_PI;
  if ( name == "e" ) return std::exp(1);
  if ( name == "I" || name == "i" ) return complex_t( 0, 1 );
  if ( mps != nullptr ) {
    auto it = mps->find(name);
    if ( it != nullptr ) return MinuitParameterLink( it );
    else if ( mps->find(name+"_Re") != nullptr && mps->find(name+"_Im") != nullptr ) {
      if( m_isCartesian ) return MinuitParameterLink( mps->find(name+"_Re") ) + 1i * MinuitParameterLink( mps->find(name+"_Im") );
      else return MinuitParameterLink( mps->find(name+"_Re") ) * fcn::exp( m_sf * 1i *MinuitParameterLink( mps->find(name+"_Im") ) ); 
    }
    else { 
      DEBUG( "Token not understood: " << name << " [map size = " << mps->size() << "]" );
    }
  }
  return Parameter( name );
}

MinuitParameterLink::MinuitParameterLink( MinuitParameter* param ) : m_parameter( param ) {}

std::string MinuitParameterLink::to_string(const ASTResolver* resolver) const
{
  auto as_expression = dynamic_cast<const MinuitExpression*>( this->m_parameter );
  if( as_expression != nullptr ) return as_expression->expression().to_string(resolver);

  if( resolver == nullptr ) return m_parameter->name();
  if( resolver->enableCompileConstants() && m_parameter != nullptr && m_parameter->flag () == Flag::CompileTimeConstant )
    return std::to_string( m_parameter->mean() ); 
  return resolver->resolvedParameter(this);

}

std::string MinuitParameterLink::name() const {
  return m_parameter->name();
}

void MinuitParameterLink::resolve( ASTResolver& resolver ) const
{
  auto as_expression = dynamic_cast<const MinuitExpression*>( this->m_parameter );
  if( as_expression != nullptr ) return as_expression->expression().resolve(resolver);
  if( m_parameter->flag() != Flag::CompileTimeConstant ) resolver.resolve(*this);
}

complex_t MinuitParameterLink::operator()() const 
{ 
  if( m_parameter == nullptr ) ERROR("Parameter does not have end-point");
  return m_parameter->mean(); 
}

const MinuitParameter& MinuitParameterLink::param() const {
  return *m_parameter ; 
}

