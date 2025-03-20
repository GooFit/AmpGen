#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Property"

#include <boost/test/unit_test.hpp>
#include "AmpGen/Property.h"

BOOST_AUTO_TEST_CASE ( constructors_right )
{
  AmpGen::Property<unsigned> param(nullptr, "test_param", 4);
  BOOST_CHECK( param == unsigned(4) );
}

BOOST_AUTO_TEST_CASE ( constructors_tuple )
{
  AmpGen::Property<std::tuple<double, std::string>> param(nullptr, "test_param",  std::make_tuple(0.1, "hello, world") );
  auto [d,s] = std::tuple<double, std::string>(param); 
  BOOST_CHECK( d == double(0.1) );
  BOOST_CHECK( s == "hello, world" );
}

