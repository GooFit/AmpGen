#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Tensor"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/Tensor.h"

using namespace AmpGen; 

BOOST_AUTO_TEST_CASE ( constructor )
{
  Tensor A({"A_11","A_12","A_13",
            "A_21","A_22","A_23"}, Tensor::dim(2,3) );

    
}

BOOST_AUTO_TEST_CASE( multiply )
{
  Tensor::Index a,b,c;
  Tensor A({1,2,3,
            4,5,6},Tensor::dim(2,3) );
  Tensor B({7,8,
           9,10,
           11,12}, Tensor::dim(3,2) );

  Tensor C = A(a,b)*B(b,c);

  BOOST_TEST( int( std::real( C[{0,0}]()) ) == 58 );
  BOOST_TEST( int( std::real( C[{0,1}]()) ) == 64 );
  BOOST_TEST( int( std::real( C[{1,0}]()) ) == 139 );
  BOOST_TEST( int( std::real( C[{1,1}]()) ) == 154 );
 
}
