#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "ArgumentPack"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/ArgumentPack.h"

DECLARE_ARGUMENT( test_argument_string, std::string);
DECLARE_ARGUMENT( test_argument_double , double );

static const test_argument_double value;

struct Mock {
  Mock( const AmpGen::ArgumentPack& args ){
    t1 = args.getArg<test_argument_double>().val;
    t2 = args.getArg<test_argument_string>("hello world").val;
  }
  template <class ...ARGS> Mock( const ARGS&... args ) : 
    Mock( AmpGen::ArgumentPack(args...) ) {}
  double      t1;
  std::string t2;
};

BOOST_AUTO_TEST_CASE ( test_string )
{
  test_argument_string f("different value");
  BOOST_CHECK( f.val == "different value");  
}

BOOST_AUTO_TEST_CASE ( test_pack )
{
  Mock test( test_argument_string("this"), test_argument_double(12) );
  BOOST_CHECK( test.t1 == 12 );
  BOOST_CHECK( test.t2 == "this");
}
