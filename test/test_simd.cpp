#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "test_simd"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/simd/utils.h"
#include "AmpGen/MetaUtils.h"

using namespace AmpGen; 
using namespace std::complex_literals; 

#if INSTRUCTION_SET != 0 

#if INSTRUCTION_SET == INSTRUCTION_SET_AVX2d
double precision = 1e-15;
namespace simd = AmpGen::AVX2d; 
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2f
float precision = 1e-7;
namespace simd = AmpGen::AVX2f; 
#elif INSTRUCTION_SET == INSTRUCTION_SET_ARM128d
double precision = 1e-15;
namespace simd = AmpGen::ARM128d; 
#endif
using scalar_t = real_v::scalar_type;

#define test_simd( avx_function, scalar_function, data) \
{ auto v = avx_function(data); \
  for(int i =0;i!=real_v::size;++i){ BOOST_TEST( utils::at(v, i) == (scalar_function(utils::at(data,i))) ); } } 

#define test_simd_complex( avx_function, scalar_function, data) \
{ auto v = avx_function(data); \
  for(int i =0;i!=real_v::size;++i){ BOOST_TEST( std::real(utils::at(v, i)) == std::real(scalar_function(utils::at(data,i))) ); \
  BOOST_TEST( std::imag(utils::at(v, i)) == std::imag(scalar_function(utils::at(data,i))) ); } } 

BOOST_AUTO_TEST_CASE( test_log, *utf::tolerance(scalar_t(50*precision) ) ) 
{
  std::vector<scalar_t> test_v = {0.3, 0.5, 10, 7.0, 5, 3, 8, 7};
  test_simd( simd::log, std::log, simd::real_v(test_v.data()) ); 
}

BOOST_AUTO_TEST_CASE( test_complex_log, *utf::tolerance(scalar_t(precision)) )
{
  std::array<scalar_t, 8> rv = {0.3, 0.5, 10, -4.0, -7.0, -.5, 0.5, 0.5}; 
  std::array<scalar_t, 8> iv = {-3.0, -4.0, 3, 1.0, 0.1, -4.0, -4.0, -4.0};
  test_simd_complex( log, std::log, simd::complex_v(real_v(rv.data()), real_v(iv.data())) ); 
}

BOOST_AUTO_TEST_CASE( test_fmod, *utf::tolerance(scalar_t(precision)) ) 
{
  std::vector<scalar_t> a = {5.1, -5.1, 5.1, -5.1, 5.1, -5.1, 5.1, -5.1};
  std::vector<scalar_t> b = {3.0, +3.0, -3.0, -3.0, 3.0, +3.0, -3.0, -3.0};

  simd::real_v av( a.data() );
  simd::real_v bv( b.data() );

  auto modv = simd::fmod(av, bv);
  
  auto mod = modv.to_array();
  for( int i = 0; i != real_v::size; ++i ){
    BOOST_TEST( mod[i] ==  fmod(a[i],b[i])  );
  }
}

BOOST_AUTO_TEST_CASE( test_double_to_int )
{ 
  std::vector<scalar_t> a = {17.4, 19.8, 12.1, 4007.3, 12.0, 14.0, -7.0, 623.};
  auto f = real_v( a.data() ).to_int();
  auto t = simd::store( f ); 
  for( int i = 0 ; i != real_v::size ;++i ) { std::cout << "test["<<i<<"] : " << t[i] << " " << a[i] << std::endl;  BOOST_TEST( t[i] == int(a[i]) ); }
}
BOOST_AUTO_TEST_CASE( test_gather )
{
                              // 0     1      2      3      4     5      6   7 
  std::vector<double> data = { 15.4, 19.7, 121.8, -15.6, M_PI, sqrt(2), 5.7, 12 };
  std::vector<scalar_t> addr = { 0.2, 5.3, 3.8, 4.1,  2.0, 3.3,7.0, 8.0 };
  auto v = simd::gather( data.data(), simd::real_v(addr.data()) ).to_array();
  for( int i = 0 ; i != real_v::size; ++i ){ BOOST_TEST( v[i] == scalar_t( data[int(addr[i])]) ); }
  
}

BOOST_AUTO_TEST_CASE( test_trig, *utf::tolerance(scalar_t(50*precision)) )
{
  std::vector<scalar_t> data = {0.1,0.4,-2.0,5.0,0.7,0.9,-53.0, 12.0};
  test_simd( simd::cos, std::cos, simd::real_v(data.data()));
  test_simd( simd::sin, std::sin, simd::real_v(data.data()));
  test_simd( simd::tan, std::tan, simd::real_v(data.data()));
}
#else 
BOOST_AUTO_TEST_CASE( test_dummy )
{
  BOOST_TEST( 1 == 1 );
}
#endif


