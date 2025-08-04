#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Property"

#include <boost/test/unit_test.hpp>
#include "AmpGen/Property.h"

BOOST_AUTO_TEST_CASE ( constructors_right )
{
  AmpGen::Property<unsigned> param(nullptr, "test_param", 4);
  BOOST_CHECK( param == unsigned(4) );
}

