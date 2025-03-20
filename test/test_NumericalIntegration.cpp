#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "NumericalIntegration"

#include "AmpGen/NumericalIntegration.h"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;
using namespace AmpGen;

BOOST_AUTO_TEST_CASE( test_NumericalIntegration )
{
  auto f = []( const std::array<real_v, 2>& v){return exp( -0.5* ( v[0]*v[0] + v[1]*v[1]  ) );};
  double v = integrate<2>( f, std::array<double, 2>{-5,-5}, std::array<double,2>{5,5} ); 
  // std::cout << v - pow( 0.5 * ( 1 + std::erf(5/sqrt(2)) ) * sqrt( 2 * M_PI ) , 2 )  << std::endl; 
  std::cout << v - 2*M_PI << std::endl; //  pow( 0.5 * ( 1 + std::erf(5/sqrt(2)) ) * sqrt( 2 * M_PI ) , 2 )  << std::endl; 
  BOOST_TEST( v ==  pow( 0.5 * ( 1 + std::erf(5/sqrt(2)) ) * sqrt( 2 * M_PI ) , 2 ), boost::test_tools::tolerance(1e-6) );
}
