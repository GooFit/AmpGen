
#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "test_avx2"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;


#if ENABLE_AVX2d
#include "AmpGen/simd/utils.h"

using namespace AmpGen; 
using namespace AmpGen::AVX2d;
using namespace std::complex_literals; 

#define test_simd( avx_function, scalar_function, data, tv) \
{ auto r = avx_function( data ).to_array(); auto vals = data.to_array(); \
  for(int i =0;i!=4;++i) BOOST_TEST( r[i] == scalar_function(vals[i]), boost::test_tools::tolerance(tv) ); } 

BOOST_AUTO_TEST_CASE( test_log )
{
  test_simd( AVX2d::log, std::log, AVX2d::real_v(0.3, 0.5, 10, 7.0), 1e-12 ); 
}

BOOST_AUTO_TEST_CASE( test_complex_log )
{
  std::array<std::complex<double>, 4> pr = {0.3 - 3.0*1i, 0.5 - 4.0*1i, 10.+3.*1i, -4.0 + 1.0*1i};
  test_simd( AVX2d::log, std::log, AVX2d::complex_v( pr.data() ), 1e-8 );
}

BOOST_AUTO_TEST_CASE( test_fmod ) 
{
  std::vector<double> a = {5.1, -5.1, 5.1, -5.1};
  std::vector<double> b = {3.0, +3.0, -3.0, -3.0};

  AVX2d::real_v av( a.data() );
  AVX2d::real_v bv( b.data() );

  auto modv = AVX2d::fmod(av,bv);
  
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
  test_simd( AVX2d::cos, std::cos, data, 1e-15);
  test_simd( AVX2d::sin, std::sin, data, 1e-15);
  test_simd( AVX2d::tan, std::tan, data, 1e-15);
}

#else 
BOOST_AUTO_TEST_CASE( test_dummy )
{
  BOOST_TEST( 1 == 1 );
}
#endif


