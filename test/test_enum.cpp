
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "enum"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/enum.h"

namespace AmpGen {
  make_enum( test_enum, state1, state2, state3 )
}
using namespace AmpGen;

BOOST_AUTO_TEST_CASE( test_enums )
{
  BOOST_CHECK( parse<test_enum> ("state1") == test_enum::state1  );
  BOOST_CHECK( parse<test_enum> ("state2") == test_enum::state2  );
  BOOST_CHECK( parse<test_enum> ("state3") == test_enum::state3  );
  BOOST_CHECK( parse<test_enum> ("blag")   == test_enum::Invalid );
}
