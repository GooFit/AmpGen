#ifndef AMPGEN_AVX_TYPES
#define AMPGEN_AVX_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#include "AmpGen/simd/avx_mathfun.h"
#include <omp.h>

namespace AmpGen {
  namespace AVX2 {
    struct real_v {
      __m256 data;
      static constexpr unsigned size = 8 ;
      typedef float scalar_type;
      real_v() = default; 
      real_v(__m256 data ) : data(data) {}
      real_v(const float& f )      : data( _mm256_set1_ps(f) ) {}
      real_v(const double& f )     : data( _mm256_set1_ps( float(f) )) {}
      real_v(const float* f )      : data( _mm256_loadu_ps( f ) ) {}
      void store( float* ptr ) const { _mm256_storeu_ps( ptr, data ); }
      std::array<float, 8> to_array() const { std::array<float, 8> b; store( &b[0] ); return b; }
      float at(const unsigned i) const { return to_array()[i] ; }       
      operator __m256() const { return data ; } 
    };
    
    inline real_v operator+( const real_v& lhs, const real_v& rhs ) { return _mm256_add_ps(lhs, rhs); }
    inline real_v operator-( const real_v& lhs, const real_v& rhs ) { return _mm256_sub_ps(lhs, rhs); }
    inline real_v operator*( const real_v& lhs, const real_v& rhs ) { return _mm256_mul_ps(lhs, rhs); }
    inline real_v operator/( const real_v& lhs, const real_v& rhs ) { return _mm256_div_ps(lhs, rhs); }
    inline real_v operator-( const real_v& x ) { return -1.f * x; }
    inline real_v operator&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
    inline real_v operator|( const real_v& lhs, const real_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
    inline real_v operator^( const real_v& lhs, const real_v& rhs ) { return _mm256_xor_ps( lhs, rhs ); }
    inline real_v operator+=(real_v& lhs, const real_v& rhs ){ lhs = lhs + rhs; return lhs; }
    inline real_v operator-=(real_v& lhs, const real_v& rhs ){ lhs = lhs - rhs; return lhs; }
    inline real_v operator*=(real_v& lhs, const real_v& rhs ){ lhs = lhs * rhs; return lhs; }
    inline real_v operator/=(real_v& lhs, const real_v& rhs ){ lhs = lhs / rhs; return lhs; }
    inline real_v operator&&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
    inline real_v operator||( const real_v& lhs, const real_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
    inline real_v operator!( const real_v& x ) { return x ^ _mm256_castsi256_ps( _mm256_set1_epi32( -1 ) ); }
    inline real_v operator<( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_LT_OS ); }
    inline real_v operator>( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_GT_OS ); }
    inline real_v operator==( const real_v& lhs, const real_v& rhs ){ return _mm256_cmp_ps( lhs, rhs, _CMP_EQ_OS ); }
    inline real_v sqrt( const real_v& v ) { return _mm256_sqrt_ps(v); } 
    inline real_v sin( const real_v& v )  { return sin256_ps(v) ; }
    inline real_v cos( const real_v& v )  { return cos256_ps(v) ; }
    inline real_v tan( const real_v& v )  { real_v s; real_v c; sincos256_ps(v, (__m256*)&s, (__m256*)&c) ; return s/c; }
    inline real_v log( const real_v& v )  { return log256_ps(v) ; }
    inline real_v exp( const real_v& v )  { return exp256_ps(v) ; }
    inline real_v abs ( const real_v& v ) { return v & _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ); }
    inline real_v select(const real_v& mask, const real_v& a, const real_v& b ) { return _mm256_blendv_ps( b, a, mask ); }
    inline real_v select(const bool& mask   , const real_v& a, const real_v& b ) { return mask ? a : b; } 
    inline real_v atan2( const real_v& y, const real_v& x ){ 
      std::array<float, 8> bx{x.to_array()}, by{y.to_array()}, rt;
      for( unsigned i = 0 ; i != 8 ; ++i ) rt[i] = std::atan2( by[i] , bx[i] );
      return real_v (rt.data() ); 
    }
    inline real_v fmadd( const real_v& a, const real_v& b, const real_v& c )
    {
      return _mm256_fmadd_ps(a, b, c );
    }
    struct complex_v {
      real_v re;
      real_v im;
      typedef std::complex<float> scalar_type;
      static constexpr unsigned size = 8 ;

      real_v real() const { return re; }
      real_v imag() const { return im; }
      real_v norm() const { return re*re + im *im ; }
      complex_v() = default; 
      complex_v( const real_v& re, const real_v& im) : re(re), im(im) {}
      complex_v( const float&   re, const float& im) : re(re), im(im) {}
      complex_v( const std::complex<double>& f ) : re( f.real() ), im( f.imag() ) {}
      complex_v( const std::complex<float>& f  ) : re( f.real() ), im( f.imag() ) {}
      const std::complex<float> at(const unsigned i) const { return std::complex<float>(re.to_array()[i], im.to_array()[i]) ; }       
      void store( float* sre, float* sim ){ re.store(sre); im.store(sim);  } 
      void store( std::complex<float>* r ){ 
        auto re_arr = re.to_array();
        auto im_arr = im.to_array();
        for( unsigned i = 0 ; i != re_arr.size(); ++i ) r[i] = std::complex<float>( re_arr[i], im_arr[i] ); 
      }
    };

    inline std::ostream& operator<<( std::ostream& os, const real_v& obj ) { 
      auto buffer = obj.to_array();
      for( unsigned i = 0 ; i != 8; ++i ) os << buffer[i] << " ";
      return os; 
    }
    inline real_v     real(const complex_v& arg ){ return arg.re ; }
    inline real_v     imag(const complex_v& arg ){ return arg.im ; }
    inline complex_v   conj(const complex_v& arg ){ return complex_v(arg.re, -arg.im) ; }
    inline real_v     conj(const real_v& arg ){ return arg ; } 
    inline complex_v operator+( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.re + rhs, lhs.im); }
    inline complex_v operator-( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.re - rhs, lhs.im); }
    inline complex_v operator*( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.re*rhs, lhs.im*rhs); }
    inline complex_v operator/( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.re/rhs, lhs.im/rhs); }
    inline complex_v operator+( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs + rhs.re,  rhs.im); }
    inline complex_v operator-( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs - rhs.re, - rhs.im); }
    inline complex_v operator*( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs*rhs.re, lhs*rhs.im); }
    inline complex_v operator/( const real_v& lhs, const complex_v& rhs ) { return complex_v( lhs * rhs.re , -lhs *rhs.im) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
    inline complex_v operator+( const complex_v& lhs, const complex_v& rhs ) { return complex_v(lhs.re + rhs.re, lhs.im + rhs.im); }
    inline complex_v operator-( const complex_v& lhs, const complex_v& rhs ) { return complex_v(lhs.re - rhs.re, lhs.im - rhs.im); }
    inline complex_v operator*( const complex_v& lhs, const complex_v& rhs ) { return complex_v(lhs.re*rhs.re - lhs.im*rhs.im, lhs.re*rhs.im  + lhs.im*rhs.re); }
    inline complex_v operator/( const complex_v& lhs, const complex_v& rhs ) { return complex_v(lhs.re*rhs.re + lhs.im*rhs.im, -lhs.re*rhs.im  + lhs.im*rhs.re) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
    inline complex_v operator-( const complex_v& x ) { return -1.f * x; }
    inline real_v abs( const complex_v& v ) { return sqrt( v.re * v.re + v.im * v.im ) ; }
    inline real_v norm( const complex_v& v ) { return  ( v.re * v.re + v.im * v.im ) ; }
    inline complex_v select(const real_v& mask, const complex_v& a, const complex_v& b ) { return complex_v( _mm256_blendv_ps( b.re, a.re, mask ), _mm256_blendv_ps( b.im, a.im, mask ) ); }
    inline complex_v select(const real_v& mask, const real_v&   a, const complex_v& b ) { return complex_v( _mm256_blendv_ps( b.re, a   , mask ), _mm256_blendv_ps( b.im, real_v(0.f), mask ) ); }
    inline complex_v select(const real_v& mask, const complex_v& a, const real_v& b   ) { return complex_v( _mm256_blendv_ps( b, a.re, mask ), _mm256_blendv_ps( real_v(0.f), a.im, mask ) ); }
    inline complex_v select(const bool& mask   , const complex_v& a, const complex_v& b ) { return mask ? a : b; }
    inline complex_v exp( const complex_v& v ){ 
      real_v s; real_v c; sincos256_ps(v.im, (__m256*)&s, (__m256*)&c) ; 
      return exp( v.re ) * complex_v(c, s);
    }
    inline std::ostream& operator<<( std::ostream& os, const complex_v& obj ) { return os << "( "<< obj.re << ") (" << obj.im << ")"; }
    #pragma omp declare reduction(+: real_v: \
	omp_out = omp_out + omp_in) 
    #pragma omp declare reduction(+: complex_v: \
	omp_out = omp_out + omp_in) 
  
  }  
}

#endif
