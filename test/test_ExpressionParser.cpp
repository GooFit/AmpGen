#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "ExpressionParser"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/ExpressionParser.h"
#include "AmpGen/MinuitParameterSet.h"

using namespace AmpGen; 

BOOST_AUTO_TEST_CASE( test_ExpressionParser ) {

  double a = 1.0;
  double b = 0.5; 

  MinuitParameterSet mps( {
      new MinuitParameter("a",MinuitParameter::Flag::Float, a,0. )
  ,   new MinuitParameter("b",MinuitParameter::Flag::Float, b,0. )
  });
  ExpressionParser::getMe()->setMPS( &mps );

  double pi = M_PI;

  auto test = [](const std::string& expr) -> double { 
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
    return std::real( ExpressionParser::Parse(newLine)() ); }; 

  BOOST_CHECK( test("a+(cos(b-sin(2/a*pi))-sin(a-cos(2*b/pi)))-b") == a+(cos(b-sin(2/a*pi))-sin(a-cos(2*b/pi)))-b );
  BOOST_CHECK( test("sin(a)+sin(b)")                               == sin(a)+sin(b) );                                  
  BOOST_CHECK( test("abs(sin(sqrt(a^2+b^2))*255)")                 == abs(sin(sqrt(a*a+b*b))*255) );
  BOOST_CHECK( test("sqrt(a)<sin(8)")                              == (sqrt(a)<sin(8)) );
  BOOST_CHECK( test("(10+sqrt(a))<(sin(8)^2)")                     == ((10+sqrt(a))<pow(sin(8),2)) );
  BOOST_CHECK( test("(b+a/b) * (a-b/a)")                           == (b+a/b) * (a-b/a) );
  BOOST_CHECK( test("(0.1*a+1)*a+1.1-sin(a)-log(a)/a*3/4")         == (0.1*a+1)*a+1.1-sin(a)-log(a)/a*3/4 );
  BOOST_CHECK( test("sin(2 * a) + cos(pi / b)")                    == sin(2 * a) + cos(pi / b) );
//  BOOST_CHECK( test("1 - sin(2 * a) + cos(pi / b)")                == 
//  BOOST_CHECK( test("sqrt(1 - sin(2 * a) + cos(pi / b) / 3)")      == 
//  BOOST_CHECK( test("(a^2 / sin(2 * pi / b)) -a / 2")              == 
//  BOOST_CHECK( test("1-(a/b*0.5)")                                 == 
//  BOOST_CHECK( test("e^log(7*a)")                                  ==
//  BOOST_CHECK( test("10^log(3+b)")                                 ==
//  BOOST_CHECK( test("(cos(2.41)/b)")                               == 
//  BOOST_CHECK( test("-(sin(pi+a)+1)")                              == 
//  BOOST_CHECK( test("a-(e^(log(7+b)))")                            == 

}
