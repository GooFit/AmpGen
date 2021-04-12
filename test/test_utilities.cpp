#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Utilities"
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/Utilities.h"

BOOST_AUTO_TEST_CASE( test_swap_char )
{
  std::string test = "hello, world"; 
  AmpGen::swapChars( test, 'h', 'w' );
  BOOST_TEST( test == "wello, horld");
}

BOOST_AUTO_TEST_CASE( test_split )
{
  std::string test = "hello world tokens";
  auto tokens = AmpGen::split(test,' ');
  BOOST_TEST(tokens[0] == "hello");
  BOOST_TEST(tokens[1] == "world");
  BOOST_TEST(tokens[2] == "tokens");
}
