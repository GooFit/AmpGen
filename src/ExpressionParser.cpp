#include "AmpGen/ExpressionParser.h"

#include <cmath>
#include <ostream>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

DEFINE_CAST( MinuitParameterLink )
DEFINE_CAST( ExpressionPack )

void ExpressionParser::processBinaryOperators( std::vector<std::string>& opCodes, std::vector<Expression>& expressions )
{
  for ( auto& fcn : m_binaryFunctions ) {
    for ( int pos = 0; pos < (int)opCodes.size(); ++pos ) {
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
  for ( int pos = 0; pos < (int)opCodes.size(); ++pos ) {
    auto fcn = m_unaryFunctions.find( opCodes[pos] );
    INFO( " op = " << opCodes[pos] << " pos = " << pos );
    if ( fcn == m_unaryFunctions.end() ) continue;
    expressions[pos] = ( fcn )->second( expressions[pos] );
    opCodes.erase( opCodes.begin() + pos );
    pos--;
  }
}

Expression ExpressionParser::parseTokens( const std::vector<std::string>& tokens )
{
  int nOpenedBrackets = 0;
  std::vector<std::vector<std::string>> expressions;
  expressions.push_back( std::vector<std::string>() );
  for ( auto iToken = tokens.begin(); iToken != tokens.end(); ++iToken ) {
    auto& token = *iToken;
    expressions.rbegin()->push_back( token );
    if ( token == "(" ) nOpenedBrackets++;
    if ( token == ")" ) nOpenedBrackets--;
    if ( nOpenedBrackets == 0 && iToken != tokens.end() - 1 ) expressions.push_back( std::vector<std::string>() );
  }
  for ( auto& eC : expressions ) {
    if ( eC.size() > 2 && *eC.begin() == "(" && *eC.rbegin() == ")" ) {
      eC.erase( eC.end() - 1 );
      eC.erase( eC.begin() );
    }
  }
  if ( expressions.size() == 1 ) return processEndPoint( expressions[0][0] );

  std::vector<std::string> opCodes;
  std::vector<Expression> parsedExpressions;
  for ( size_t pos = 0; pos < expressions.size(); pos++ ) {
    if ( pos % 2 == expressions.size() % 2 )
      opCodes.push_back( expressions[pos][0] );
    else
      parsedExpressions.push_back( parseTokens( expressions[pos] ) );
  }
  processUnaryOperators( opCodes, parsedExpressions );
  processBinaryOperators( opCodes, parsedExpressions );
  return parsedExpressions[0];
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

  add_binary( "/", []( auto& A, auto& B ) { return A / B; } ); ///  operator ordering here matters!
  add_binary( "*", []( auto& A, auto& B ) { return A * B; } );
  add_binary( "+", []( auto& A, auto& B ) { return A + B; } );
  add_binary( "-", []( auto& A, auto& B ) { return A - B; } );
  add_binary( ">", []( auto& A, auto& B ) { return A > B; } );
  add_binary( "<", []( auto& A, auto& B ) { return A < B; } );
  add_binary( "&&", []( auto& A, auto& B ) { return A && B; } );
  add_binary( ",", []( auto& A, auto& B ) { return ExpressionPack( A, B ); } );
}

Expression ExpressionParser::Parse( const std::string& str ) { return getMe()->parseTokens( split( str, ' ' ) ); }
ExpressionParser* ExpressionParser::gExpressionParser = nullptr;

Expression ExpressionParser::processEndPoint( const std::string& name )
{
  bool status  = true;
  double value = lexical_cast<double>( name, status );
  if ( status == true ) return Constant( value );
  if ( name == "PI" ) return Constant( M_PI );
  if ( name == "I" ) return Constant( 0, 1 );

  if ( m_mps != nullptr ) {
    auto map = m_mps->map();
    auto it  = map.find( name );
    if ( it != map.end() )
      return MinuitParameterLink( it->second );
    else {
      WARNING( "Token not understood: " << name << " [map size = " << map.size() << "]" );
      for ( auto& ip : map ) INFO( "map entry = " << ip.first );
    }
  }
  return Parameter( name, 0, true );
}
MinuitParameterLink::MinuitParameterLink( MinuitParameter* param ) : m_parameter( param ) {}
std::string MinuitParameterLink::to_string() const { return m_parameter->name(); }

void MinuitParameterLink::resolve( ASTResolver& resolver ){}
std::complex<double> MinuitParameterLink::operator()() const { return m_parameter->mean(); }

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
std::string ExpressionPack::to_string() const
{
  std::string rt = "";
  for ( auto expr : m_expressions ) {
    rt += expr.to_string() + ", ";
  }
  return rt.substr( 0, rt.length() - 2 );
}

void ExpressionPack::resolve( ASTResolver& resolver )
{
  for ( auto& expr : m_expressions ) expr.resolve( resolver );
}
std::complex<double> ExpressionPack::operator()() const { return 0; }
