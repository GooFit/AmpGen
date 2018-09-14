#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "NamedParameter"

#include <boost/test/unit_test.hpp>
#include "AmpGen/NamedParameter.h"

BOOST_AUTO_TEST_CASE ( constructors_right )
{
  AmpGen::NamedParameter<size_t> param("test_param", 4);
  BOOST_CHECK( param == 4 );
}

