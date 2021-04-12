#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Expression"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/Expression.h"

BOOST_AUTO_TEST_CASE ( sum_to_string )
{
  AmpGen::Parameter A("A");
  AmpGen::Parameter B("B");
  BOOST_CHECK( (A+B).to_string() == "A + B" );
}

BOOST_AUTO_TEST_CASE ( sum_evaluate )
{
  AmpGen::Parameter A("A",4);
  AmpGen::Parameter B("B",5);
  BOOST_CHECK( (A+B)() == std::complex<double>(9,0) );
}

BOOST_AUTO_TEST_CASE ( product_to_string )
{
  AmpGen::Parameter A("A");
  AmpGen::Parameter B("B");
  BOOST_CHECK( (A*B).to_string() == "A*B" );
}

BOOST_AUTO_TEST_CASE ( product_evaluate )
{
  AmpGen::Parameter A("A",4);
  AmpGen::Parameter B("B",5);
  BOOST_CHECK( (A*B)() == std::complex<double>(20,0) );
}

BOOST_AUTO_TEST_CASE( test_composite, * utf::tolerance(1e-6) )
{
  AmpGen::Parameter A("A",4);
  AmpGen::Parameter B("B",5);
  double value = std::real( ( A*B/( AmpGen::fcn::cos(B) * AmpGen::fcn::sqrt(A) ))() );
  BOOST_TEST( value == (4*5/(cos(5)*sqrt(4))), boost::test_tools::tolerance(1e-6)) ;
}
