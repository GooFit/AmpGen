#include "AmpGen/ExpressionParser.h"

#include <stddef.h>
#include <cmath>
#include <ostream>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/ASTResolver.h"

using namespace AmpGen;
using namespace std::complex_literals; 

DEFINE_CAST( MinuitParameterLink )
DEFINE_CAST( ExpressionPack )

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
    expressions[pos] = ( fcn )->second( expressions[pos] );
    opCodes.erase( opCodes.begin() + pos );
    pos--;
  }
}

std::vector<std::string>::iterator findMatchingBracket( 
   const std::vector<std::string>::iterator& begin, 
   const std::vector<std::string>::iterator& end )
{
  int openedBrackets = 1;
  if( begin + 1 >= end ) return end; 
  for( auto it = begin+1; it != end; ++it )
  {
    if( *it == "(")  openedBrackets++;
    if( *it == ")" ) openedBrackets--;
    if( openedBrackets == 0 ) return it;
  }
  ERROR("Unmatched brace in expression");
  return end; 
}

Expression ExpressionParser::parseTokens( const std::vector<std::string>::iterator& begin,
                                          const std::vector<std::string>::iterator& end )
{
  std::vector<std::string> opCodes;
  std::vector<Expression> expressions; 
  for( auto it = begin; it != end; ++it )
  {
    if( *it == "(" ){
      auto begin_2 = it;
      auto end_2   = findMatchingBracket(it, end);
      expressions.emplace_back( parseTokens(begin_2+1, end_2) );
      it = end_2;
    }
    else {
      auto f = std::find_if(m_binaryFunctions.begin(), m_binaryFunctions.end(), [it](auto& jt){ return jt.first == *it; } );
      if( f != m_binaryFunctions.end() || m_unaryFunctions.count(*it) ) opCodes.push_back(*it);
      else expressions.push_back( processEndPoint( *it ) );
    }
  }
  processUnaryOperators( opCodes, expressions );
  processBinaryOperators( opCodes, expressions );
  if( expressions.size() != 1 ){
    ERROR("Could not process expression!");
  }
  return expressions[0];
}

ExpressionParser::ExpressionParser() : m_mps( nullptr )
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

  add_binary( "^" , [](auto& A, auto& B ) { return fcn::pow(A,B); } );
  add_binary( "/" , [](auto& A, auto& B ) { return A / B; } ); 
  add_binary( "*" , [](auto& A, auto& B ) { return A * B; } );
  add_binary( "+" , [](auto& A, auto& B ) { return A + B; } );
  add_binary( "-" , [](auto& A, auto& B ) { return A - B; } );
  add_binary( ">" , [](auto& A, auto& B ) { return A > B; } );
  add_binary( "<" , [](auto& A, auto& B ) { return A < B; } );
  add_binary( "&&", []( auto& A, auto& B ) { return A && B; } );
  add_binary( "," , []( auto& A, auto& B ) { return ExpressionPack( A, B ); } );
}

Expression ExpressionParser::Parse( const std::string& str ) { 
  auto tokens = split(str, ' ');
  return getMe()->parseTokens(tokens.begin(), tokens.end() ); 
}
ExpressionParser* ExpressionParser::gExpressionParser = nullptr;

Expression ExpressionParser::processEndPoint( const std::string& name )
{
  bool status  = true;
  double value = lexical_cast<double>( name, status );
  if ( status == true ) return value;
  if ( name == "PI" ) return M_PI;
  if ( name == "pi" ) return M_PI;
  if ( name == "e" ) return std::exp(1);
  if ( name == "I" ) return complex_t( 0, 1 );
  if ( name == "i" ) return complex_t( 0, 1 );
  if ( m_mps != nullptr ) {
    auto it = m_mps->find(name);
    if ( it != nullptr ) return MinuitParameterLink( it );
    else if ( m_mps->find(name+"_Re") != nullptr && m_mps->find(name+"_Im") != nullptr ) {
      return MinuitParameterLink( m_mps->find(name+"_Re") ) + 1i * MinuitParameterLink( m_mps->find(name+"_Im") );
    }
    else { 
      WARNING( "Token not understood: " << name << " [map size = " << m_mps->size() << "]" );
      for( auto& m : *m_mps ) INFO( m->name() );
    }
  }
  return Parameter( name, 0, true );
}

MinuitParameterLink::MinuitParameterLink( MinuitParameter* param ) : m_parameter( param ) {}
std::string MinuitParameterLink::to_string(const ASTResolver* resolver) const
{
  return resolver == nullptr ? m_parameter->name() : resolver->resolvedParameter(this);
}

std::string MinuitParameterLink::name() const {
  return m_parameter->name();
}

void MinuitParameterLink::resolve( ASTResolver& resolver ) const
{
  resolver.resolve(*this);
}

complex_t MinuitParameterLink::operator()() const { 
  if( m_parameter == nullptr ) ERROR("Parameter does not have end-point");
  return m_parameter->mean(); }

ExpressionPack::ExpressionPack( const Expression& A, const Expression& B )
{
  auto c1 = dynamic_cast<ExpressionPack*>( A.get() );
  auto c2 = dynamic_cast<ExpressionPack*>( B.get() );
  if ( c1 != nullptr ) {
    for ( auto& expr : c1->m_expressions ) m_expressions.push_back( expr );
  } else
    m_expressions.push_back( A );
  if ( c2 != nullptr ) {
    for ( auto& expr : c2->m_expressions ) m_expressions.push_back( expr );
  } else
    m_expressions.push_back( B );
}
std::string ExpressionPack::to_string(const ASTResolver* resolver) const
{
  std::string rt = "";
  for ( auto expr : m_expressions ) {
    rt += expr.to_string(resolver) + ", ";
  }
  return rt.substr( 0, rt.length() - 2 );
}

void ExpressionPack::resolve( ASTResolver& resolver ) const
{
  for ( auto& expr : m_expressions ) expr.resolve( resolver );
}

complex_t ExpressionPack::operator()() const { return 0; }
