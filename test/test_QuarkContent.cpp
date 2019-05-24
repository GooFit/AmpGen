
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "QuarkContent"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/QuarkContent.h"

using namespace AmpGen; 

BOOST_AUTO_TEST_CASE( test_QuarkState ) {

  QuarkState uD("uD");
  QuarkState dU("dU");
  
  BOOST_CHECK( uD != dU );

  auto vu = uD + dU;

  BOOST_CHECK( vu.isVacuum());

}
