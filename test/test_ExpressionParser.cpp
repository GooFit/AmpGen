#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "ExpressionParser"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/ExpressionParser.h"
#include "AmpGen/MinuitParameterSet.h"

using namespace AmpGen; 
  
BOOST_AUTO_TEST_CASE( simple_numericalExpressions ){
  auto test = [](const std::string& expr) -> double { return std::real( ExpressionParser::parse(expr)() ); }; 
  BOOST_CHECK( test("1 + 2 * 3"                             ) ==  1 + 2 * 3 );
  BOOST_CHECK( test("4 + 2                          "       ) ==  6 ); 
  BOOST_CHECK( test("3 * 6 - 7 + 2                  "       ) ==  3 * 6 - 7 + 2 );
  BOOST_CHECK( test("6 * 2 + ( 5 - 3 ) * 3 - 8        "     ) ==  6 * 2 + ( 5 - 3)*3 -8 );
  BOOST_CHECK( test("( 3 + 4 ) + 7 * 2 - 1 - 9        "     ) ==  3 + 4 + 7 * 2 - 1 - 9);
  BOOST_CHECK( test("5 - 2 + 4 * ( 8 - ( 5 + 1 ) ) + 9  "   ) ==  5 - 2 + 4 *( 8 - (5+1)) +9 );
  BOOST_CHECK( test("( 8 - 1 + 3 ) * 6 - ( ( 3 + 7 ) * 2 )" ) == (8-1+3)*6 - ((3+7)*2) );
}

BOOST_AUTO_TEST_CASE( parametericExpressions ) {

  double a = 1.0;
  double b = 0.5; 
  double c = 4.0; 
  double d = 0.2;
  double f = 3.13;
  
  MinuitParameterSet mps( {
    new MinuitParameter("a", Flag::Free, a, 0.)
  , new MinuitParameter("b", Flag::Free, b, 0.)
  , new MinuitParameter("c", Flag::Free, c, 0.)
  , new MinuitParameter("d", Flag::Free, d, 0.)
  , new MinuitParameter("f", Flag::Free, f, 0.)
  });

  double pi = M_PI;

  auto test = [mps](const std::string& expr) -> double { 
    std::string newLine=""; 
    for( auto& ch : expr ){
    if( ch == '(' ||
        ch == ')' ||
        ch == '+' ||
        ch == '-' ||
        ch == '*' ||
        ch == '/' ||
        ch == '>' ||
        ch == '<' ||
        ch == '^' ){
      newLine.push_back(' ');
      newLine.push_back(ch);
      newLine.push_back(' ');
    }
    else newLine.push_back(ch);
  }
  return std::real( ExpressionParser::parse(newLine, &mps)() ); }; 

  BOOST_CHECK( test("a+(cos(b-sin(2/a*pi))-sin(a-cos(2*b/pi)))-b") == a+(cos(b-sin(2/a*pi))-sin(a-cos(2*b/pi)))-b );
  BOOST_CHECK( test("sin(a)+sin(b)")                               == sin(a)+sin(b) );                                  
  BOOST_CHECK( test("abs(sin(sqrt(a^2+b^2))*255)")                 == abs(sin(sqrt(a*a+b*b))*255) );
  BOOST_CHECK( test("sqrt(a)<sin(8)")                              == (sqrt(a)<sin(8)) );
  BOOST_CHECK( test("(10+sqrt(a))<(sin(8)^2)")                     == ((10+sqrt(a))<pow(sin(8),2)) );
  BOOST_CHECK( test("(b+a/b) * (a-b/a)")                           == (b+a/b) * (a-b/a) );
  BOOST_CHECK( test("(0.1*a+1)*a+1.1-sin(a)-log(a)/a*3/4")         == (0.1*a+1)*a+1.1-sin(a)-log(a)/a*3/4 );
  BOOST_CHECK( test("sin(2 * a) + cos(pi / b)")                    == sin(2 * a) + cos(pi / b) );
  BOOST_CHECK( test("a - ( b + c ) "                        )      ==  a - ( b + c ) );
  BOOST_CHECK( test("a - b + c - d / b + f - a "            )      ==  a - b + c - d / b + f -a  );
}


