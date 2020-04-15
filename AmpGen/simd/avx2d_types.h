#ifndef AMPGEN_AVXd_TYPES
#define AMPGEN_AVXd_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#include "AmpGen/simd/avx_mathfun.h"
#include <omp.h>

namespace AmpGen {
  namespace AVX2d {
    #define stl_fallback( x ) \
      inline float_t x( const float_t& v ){  auto a = v.to_array(); return float_t( std::x(a[0]), std::x(a[1]), std::x(a[2]), std::x(a[3]) ) ; } 
       
    struct float_t {
      __m256d data;
      static constexpr unsigned size = 4;
      typedef double scalar_type;
      float_t() = default; 
      float_t(__m256d data ) : data(data) {}
      float_t(const double& f ) : data( _mm256_set1_pd( f )) {}
      float_t(const double& x0, const double& x1, const double& x2, const double& x3 )
      {
        double tmp[4] = {x0,x1,x2,x3};
        _mm256_loadu_pd(tmp); 
      }
      float_t(const double* f ) : data( _mm256_loadu_pd( f ) ) {}
      void store( double* ptr ) const { _mm256_storeu_pd( ptr, data ); }
      std::array<double, 4> to_array() const { std::array<double, 4> b; store( &b[0] ); return b; }
      double at(const unsigned i) const { return to_array()[i] ; }       
      operator __m256d() const { return data ; } 
    };

    inline float_t operator+( const float_t& lhs, const float_t& rhs ) { return _mm256_add_pd(lhs, rhs); }
    inline float_t operator-( const float_t& lhs, const float_t& rhs ) { return _mm256_sub_pd(lhs, rhs); }
    inline float_t operator*( const float_t& lhs, const float_t& rhs ) { return _mm256_mul_pd(lhs, rhs); }
    inline float_t operator/( const float_t& lhs, const float_t& rhs ) { return _mm256_div_pd(lhs, rhs); }
    inline float_t operator-( const float_t& x ) { return -1.f * x; }
    inline float_t operator&( const float_t& lhs, const float_t& rhs ) { return _mm256_and_pd( lhs, rhs ); }
    inline float_t operator|( const float_t& lhs, const float_t& rhs ) { return _mm256_or_pd( lhs, rhs ); }
    inline float_t operator^( const float_t& lhs, const float_t& rhs ) { return _mm256_xor_pd( lhs, rhs ); }
    inline float_t operator+=(float_t& lhs, const float_t& rhs ){ lhs = lhs + rhs; return lhs; }
    inline float_t operator-=(float_t& lhs, const float_t& rhs ){ lhs = lhs - rhs; return lhs; }
    inline float_t operator*=(float_t& lhs, const float_t& rhs ){ lhs = lhs * rhs; return lhs; }
    inline float_t operator/=(float_t& lhs, const float_t& rhs ){ lhs = lhs / rhs; return lhs; }
    inline float_t operator&&( const float_t& lhs, const float_t& rhs ) { return _mm256_and_pd( lhs, rhs ); }
    inline float_t operator||( const float_t& lhs, const float_t& rhs ) { return _mm256_or_pd( lhs, rhs ); }
    inline float_t operator!( const float_t& x ) { return x ^ _mm256_castsi256_pd( _mm256_set1_epi32( -1 ) ); }
    inline float_t operator<( const float_t& lhs, const float_t& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_LT_OS ); }
    inline float_t operator>( const float_t& lhs, const float_t& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_GT_OS ); }
    inline float_t operator==( const float_t& lhs, const float_t& rhs ){ return _mm256_cmp_pd( lhs, rhs, _CMP_EQ_OS ); }
    inline float_t sqrt( const float_t& v ) { return _mm256_sqrt_pd(v); } 
    inline float_t abs ( const float_t& v ) { return _mm256_andnot_pd(_mm256_set1_pd(-0.), v);  }
    // inline float_t sin( const float_t& v )  { return sin256_pd(v) ; }
    // inline float_t cos( const float_t& v )  { return cos256_pd(v) ; }
    // inline float_t tan( const float_t& v )  { float_t s; float_t c; sincos256_pd(v, (__m256*)&s, (__m256*)&c) ; return s/c; }
    // inline float_t log( const float_t& v )  { return log256_ps(v) ; }
    // inline float_t exp( const float_t& v )  { return exp256_ps(v) ; }
    inline float_t select(const float_t& mask, const float_t& a, const float_t& b ) { return _mm256_blendv_pd( b, a, mask ); }
    inline float_t select(const bool& mask   , const float_t& a, const float_t& b ) { return mask ? a : b; } 
    inline float_t atan2( const float_t& y, const float_t& x ){ 
      std::array<double, 4> bx{x.to_array()}, by{y.to_array()};
      return float_t ( 
          std::atan2(   by[0], bx[0]) 
          , std::atan2( by[1], bx[1]) 
          , std::atan2( by[2], bx[2]) 
          , std::atan2( by[3], bx[3]) ); 
    }
    inline __m256i double_to_int( const float_t& x )
    {
      // based on:  https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
      return _mm256_sub_epi64(_mm256_castpd_si256(x + _mm256_set1_pd(0x0018000000000000)), 
             _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000)));
    }
    inline float_t gather( const double* base_addr, const float_t& offsets) 
    {
     return _mm256_i64gather_pd(base_addr, double_to_int(offsets),sizeof(double));
    } 
    stl_fallback( log )   
    stl_fallback( exp )
    stl_fallback( tan )
    stl_fallback( sin )
    stl_fallback( cos )
    
    inline float_t remainder( const float_t& a, const float_t& b ){ return a - _mm256_round_pd(a/b, _MM_FROUND_TO_NEG_INF) * b; }
    inline float_t fmod( const float_t& a, const float_t& b )
    {
      auto r = remainder( abs(a), abs(b) );
      return select( a > 0., r, -r );
    }

    struct complex_t {
      float_t re;
      float_t im;
      typedef std::complex<double> scalar_type;

      float_t real() const { return re; }
      float_t imag() const { return im; }
      float_t norm() const { return re*re + im *im ; }
      complex_t() = default; 
      complex_t( const float_t& re, const float_t& im) : re(re), im(im) {}
      complex_t( const float&   re, const float& im) : re(re), im(im) {}
      complex_t( const std::complex<double>& f ) : re( f.real() ), im( f.imag() ) {}
      complex_t( const std::complex<float>& f  ) : re( f.real() ), im( f.imag() ) {}
      const std::complex<double> at(const unsigned i) const { return std::complex<float>(re.to_array()[i], im.to_array()[i]) ; }       
      void store( double* sre, double* sim ){ re.store(sre); im.store(sim);  } 
      void store( scalar_type* r ) const { 
        auto re_arr = re.to_array();
        auto im_arr = im.to_array();
        for( unsigned i = 0 ; i != re_arr.size(); ++i ) r[i] = scalar_type( re_arr[i], im_arr[i] ); 
      }
      auto to_array() const 
      {
        std::array<scalar_type, 4> rt;
        store( rt.data() );
        return rt;
      }
    };
   
    inline std::ostream& operator<<( std::ostream& os, const float_t& obj ) { 
      auto buffer = obj.to_array();
      for( unsigned i = 0 ; i != 4; ++i ) os << buffer[i] << " ";
      return os; 
    }
    inline float_t     real(const complex_t& arg ){ return arg.re ; }
    inline float_t     imag(const complex_t& arg ){ return arg.im ; }
    inline complex_t   conj(const complex_t& arg ){ return complex_t(arg.re, -arg.im) ; }
    inline float_t     conj(const float_t& arg ){ return arg ; } 
    inline complex_t operator+( const complex_t& lhs, const float_t& rhs ) { return complex_t(lhs.re + rhs, lhs.im); }
    inline complex_t operator-( const complex_t& lhs, const float_t& rhs ) { return complex_t(lhs.re - rhs, lhs.im); }
    inline complex_t operator*( const complex_t& lhs, const float_t& rhs ) { return complex_t(lhs.re*rhs, lhs.im*rhs); }
    inline complex_t operator/( const complex_t& lhs, const float_t& rhs ) { return complex_t(lhs.re/rhs, lhs.im/rhs); }
    inline complex_t operator+( const float_t& lhs, const complex_t& rhs ) { return complex_t(lhs + rhs.re,  rhs.im); }
    inline complex_t operator-( const float_t& lhs, const complex_t& rhs ) { return complex_t(lhs - rhs.re, - rhs.im); }
    inline complex_t operator*( const float_t& lhs, const complex_t& rhs ) { return complex_t(lhs*rhs.re, lhs*rhs.im); }
    inline complex_t operator/( const float_t& lhs, const complex_t& rhs ) { return complex_t( lhs * rhs.re , -lhs *rhs.im) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
    inline complex_t operator+( const complex_t& lhs, const complex_t& rhs ) { return complex_t(lhs.re + rhs.re, lhs.im + rhs.im); }
    inline complex_t operator-( const complex_t& lhs, const complex_t& rhs ) { return complex_t(lhs.re - rhs.re, lhs.im - rhs.im); }
    inline complex_t operator*( const complex_t& lhs, const complex_t& rhs ) { return complex_t(lhs.re*rhs.re - lhs.im*rhs.im, lhs.re*rhs.im  + lhs.im*rhs.re); }
    inline complex_t operator/( const complex_t& lhs, const complex_t& rhs ) { return complex_t(lhs.re*rhs.re + lhs.im*rhs.im, -lhs.re*rhs.im  + lhs.im*rhs.re) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
    inline complex_t operator-( const complex_t& x ) { return -1.f * x; }
    inline float_t abs( const complex_t& v ) { return sqrt( v.re * v.re + v.im * v.im ) ; }
    inline float_t norm( const complex_t& v ) { return  ( v.re * v.re + v.im * v.im ) ; }
    inline complex_t select(const float_t& mask, const complex_t& a, const complex_t& b ) { return complex_t( select(mask, a.re, b.re), select(mask, a.im, b.im ) ) ; }
    inline complex_t select(const float_t& mask, const float_t&   a, const complex_t& b ) { return complex_t( select(mask, a   , b.re), select(mask, 0.f, b.im) ); }
    inline complex_t select(const float_t& mask, const complex_t& a, const float_t& b   ) { return complex_t( select(mask, a.re, b )  , select(mask, a.im, 0.f) ); }
    inline complex_t select(const bool& mask   , const complex_t& a, const complex_t& b ) { return mask ? a : b; }
    inline complex_t exp( const complex_t& v ){ 
      return exp( v.re) * complex_t( cos( v.im ), sin( v.im ) );
    }
    inline float_t fmadd( const float_t& a, const float_t& b, const float_t& c )
    {
      return _mm256_fmadd_pd(a, b, c );
    }
    inline complex_t sqrt( const complex_t& v )
    {
      auto r = abs(v);
      return complex_t ( sqrt( 0.5 * (r + v.re) ), sqrt( 0.5*( r - v.re ) ) );
    }

    inline std::ostream& operator<<( std::ostream& os, const complex_t& obj ) { return os << "( "<< obj.re << ") (" << obj.im << ")"; }
    #pragma omp declare reduction(+: float_t: \
	omp_out = omp_out + omp_in) 
    #pragma omp declare reduction(+: complex_t: \
	omp_out = omp_out + omp_in) 
  
  }  
}

#endif
