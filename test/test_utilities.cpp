#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "ArgumentPack"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/Utilities.h"


BOOST_AUTO_TEST_CASE( test_swap_char )
{
  std::string test = "hello, world"; 
  AmpGen::swapChars( test, 'h', 'w' );
  BOOST_CHECK( test == "wello, horld");
}

