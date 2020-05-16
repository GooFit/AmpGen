
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "test_avx2"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;


#if ENABLE_AVX
#include "AmpGen/simd/utils.h"

using namespace AmpGen; 
using namespace AmpGen::AVX2d;

BOOST_AUTO_TEST_CASE( test_log )
{
  AVX2d::real_v p(0.3, 0.5, 10.0, 7.0);
  auto logged = AVX2d::log( p ).to_array() ;
  BOOST_TEST( logged[0] == std::log(0.3), boost::test_tools::tolerance(1e-12 ) );
  BOOST_TEST( logged[1] == std::log(0.5), boost::test_tools::tolerance(1e-12 ) );
  BOOST_TEST( logged[2] == std::log(10.0), boost::test_tools::tolerance(1e-12 ) );
  BOOST_TEST( logged[3] == std::log(7.0), boost::test_tools::tolerance(1e-12 ) );
}

BOOST_AUTO_TEST_CASE( test_fmod ) 
{
  std::vector<double> a = {5.1, -5.1, 5.1, -5.1};
  std::vector<double> b = {3.0, +3.0, -3.0, -3.0};

  AVX2d::real_v av( a.data() );
  AVX2d::real_v bv( b.data() );

  auto modv = AVX2d::fmod(av,bv);
  BOOST_TEST_MESSAGE( "fmod = " << modv );
  
  auto mod = modv.to_array();
  BOOST_TEST( mod[0] == 2.1  , boost::test_tools::tolerance(1e-15));
  BOOST_TEST( mod[1] == -2.1 , boost::test_tools::tolerance(1e-15));
  BOOST_TEST( mod[2] == 2.1  , boost::test_tools::tolerance(1e-15));
  BOOST_TEST( mod[3] == -2.1 , boost::test_tools::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE( test_double_to_int )
{ 
  std::vector<double> a = {17.4, 19.2, 12.1, 4007.3};
  auto f = udouble_to_uint( real_v( a.data() ));
  alignas(32) uint64_t t[ utils::size<real_v>::value ];
  _mm256_store_si256( (__m256i*)t, f);
  BOOST_TEST( t[0] == 17 );
  BOOST_TEST( t[1] == 19 );
  BOOST_TEST( t[2] == 12 );
  BOOST_TEST( t[3] == 4007 );
}

BOOST_AUTO_TEST_CASE( test_gather )
{
                              // 0     1      2      3      4     5      6 
  std::vector<double> data = { 15.4, 19.7, 121.8, -15.6, M_PI, sqrt(2), 5.7, 12 };
  std::vector<double> addr = { 0.2, 5.3, 3.1, 4.1 };
  auto v = AVX2d::gather( data.data(), AVX2d::real_v(addr.data()) ).to_array();
  BOOST_TEST( v[0] == data[0] );
  BOOST_TEST( v[1] == data[5] );
  BOOST_TEST( v[2] == data[3] );
  BOOST_TEST( v[3] == data[4] );
}

BOOST_AUTO_TEST_CASE( test_trig )
{
  auto data = AVX2d::real_v(0.1,0.4,-2.0,5.0);
  auto cos = AVX2d::cos(data).to_array();
  BOOST_TEST( cos[0] == std::cos( data.at(0 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( cos[1] == std::cos( data.at(1 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( cos[2] == std::cos( data.at(2))  , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( cos[3] == std::cos( data.at(3 )) , boost::test_tools::tolerance(1e-15) );
  
  auto sin = AVX2d::sin(data).to_array();
  BOOST_TEST( sin[0] == std::sin( data.at(0 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( sin[1] == std::sin( data.at(1 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( sin[2] == std::sin( data.at(2))  , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( sin[3] == std::sin( data.at(3 )) , boost::test_tools::tolerance(1e-15) );
  
  auto tan = AVX2d::tan(data).to_array();

  BOOST_TEST( tan[0] == std::tan( data.at(0 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( tan[1] == std::tan( data.at(1 )) , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( tan[2] == std::tan( data.at(2))  , boost::test_tools::tolerance(1e-15) );
  BOOST_TEST( tan[3] == std::tan( data.at(3 )) , boost::test_tools::tolerance(1e-15) );
}


#else 
BOOST_AUTO_TEST_CASE( test_dummy )
{
  BOOST_TEST( 1 == 1 );
}
#endif


