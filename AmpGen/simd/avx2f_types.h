#ifndef AMPGEN_AVX_TYPES
#define AMPGEN_AVX_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#include <omp.h>
#include <cmath>
#include "AmpGen/Complex.h"

#if USE_MVEC 
  extern "C" void    _ZGVdN8vvv_sincos(__m256 x, __m256i ptrs, __m256i ptrc);
#define libmvec_alias(F, O) \
  extern "C" __m256 _ZGVcN8v_##F(__m256 x);             \
  inline real_v O( const real_v& v ){ return _ZGVcN8v_##F(v) ; }
#else
#define libmvec_alias(F, O) \
  inline real_v O( const real_v& v ){ auto arr = v.to_ptr(); return real_v(  \
      std::F(arr[0]), std::F(arr[1]), std::F(arr[2]), std::F(arr[3]),          \
      std::F(arr[4]), std::F(arr[5]), std::F(arr[6]), std::F(arr[7]) ) ; }
#endif

namespace AmpGen {
  namespace AVX2f {
    struct real_v {
      __m256 data;
      static constexpr unsigned size = 8 ;
      typedef float scalar_type;
      real_v() = default; 
      real_v(__m256 data ) : data(data) {}
      real_v(const int& f ) : real_v(scalar_type(f)) {}
      real_v(const scalar_type& f ) : data( _mm256_set1_ps(f) ) {}
      real_v(const double& f )      : data( _mm256_set1_ps( scalar_type(f) )) {}
      explicit real_v(const scalar_type* f )      : data( _mm256_loadu_ps( f ) ) {}
      real_v(const scalar_type& x0, const scalar_type& x1, const scalar_type& x2, const scalar_type& x3,
          const scalar_type& x4, const scalar_type& x5, const scalar_type& x6, const scalar_type& x7)
      {
        data = _mm256_set_ps(x7,x6,x5,x4,x3,x2,x1,x0); 
      }

      void store( scalar_type* ptr ) const { _mm256_storeu_ps( ptr, data ); }
      std::array<scalar_type, 8> to_array() const { std::array<scalar_type, 8> b; store( &b[0] ); return b; }
      const scalar_type* to_ptr() const { return reinterpret_cast<const scalar_type*>( &data ) ; }
            scalar_type* to_ptr() { return reinterpret_cast<scalar_type*>( &data ) ; }
      scalar_type at(const unsigned i) const { return to_ptr()[i] ; }       
      operator __m256() const { return data ; } 
      inline real_v operator+=(const real_v& rhs ); 
      inline real_v operator-=(const real_v& rhs ); 
      inline real_v operator*=(const real_v& rhs ); 
      inline real_v operator/=(const real_v& rhs ); 
      inline __m256i to_int() const { return _mm256_cvtps_epi32(data); }
    };
    
    inline real_v operator+( const real_v& lhs, const real_v& rhs ) { return _mm256_add_ps(lhs, rhs); }
    inline real_v operator-( const real_v& lhs, const real_v& rhs ) { return _mm256_sub_ps(lhs, rhs); }
    inline real_v operator*( const real_v& lhs, const real_v& rhs ) { return _mm256_mul_ps(lhs, rhs); }
    inline real_v operator/( const real_v& lhs, const real_v& rhs ) { return _mm256_div_ps(lhs, rhs); }
    inline real_v operator-( const real_v& x ) { return -1.f * x; }
    inline real_v operator&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
    inline real_v operator|( const real_v& lhs, const real_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
    inline real_v operator^( const real_v& lhs, const real_v& rhs ) { return _mm256_xor_ps( lhs, rhs ); }
    inline real_v operator&&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
    inline real_v operator||( const real_v& lhs, const real_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
    inline real_v operator!( const real_v& x ) { return x ^ _mm256_castsi256_ps( _mm256_set1_epi32( -1 ) ); }
    inline real_v operator<( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_LT_OS ); }
    inline real_v operator>( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_GT_OS ); }
    inline real_v operator<=( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_LE_OS ); }
    inline real_v operator>=( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_GE_OS ); }
    inline real_v operator==( const real_v& lhs, const real_v& rhs ){ return _mm256_cmp_ps( lhs, rhs, _CMP_EQ_OS ); }
    inline real_v sqrt( const real_v& v ) { return _mm256_sqrt_ps(v); } 
    inline real_v real_v::operator+=(const real_v& rhs ){ *this = *this + rhs; return *this; }
    inline real_v real_v::operator-=(const real_v& rhs ){ *this = *this - rhs; return *this; }
    inline real_v real_v::operator*=(const real_v& rhs ){ *this = *this * rhs; return *this; }
    inline real_v real_v::operator/=(const real_v& rhs ){ *this = *this / rhs; return *this; }
    libmvec_alias(sinf, sin)
    libmvec_alias(cosf, cos)
    libmvec_alias(expf, exp)
    libmvec_alias(logf, log)
    inline std::array<int32_t, real_v::size> store( const __m256i& v )
    {
      alignas(32) std::array<int32_t, real_v::size> rt; 
       _mm256_store_si256( (__m256i*)&rt[0], v);
       return rt; 
    }
    
    inline void sincos( const real_v& v, real_v& s, real_v& c )
    { /// todo - add back support for mvec
      s = sin(v);
      c = cos(v);
    }
    inline std::pair<real_v, real_v> sincos( const real_v& v )
    {
      std::pair<real_v, real_v> rt;
      sincos( v, rt.first, rt.second );
      return rt; 
    }
    inline real_v tan( const real_v& v )
    {
      auto [s,c] = sincos( v );
      return s / c ;  
    }
    
    inline real_v abs   ( const real_v& v ) { return v & _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ); }
    inline real_v select(const real_v& mask, const real_v& a, const real_v& b ) { return _mm256_blendv_ps( b, a, mask ); }
    inline real_v select(const bool& mask   , const real_v& a, const real_v& b ) { return mask ? a : b; } 
    inline real_v sign  ( const real_v& v){ return select( v > 0., +1., -1. ); }
    inline real_v fmadd ( const real_v& a, const real_v& b, const real_v& c ) { return _mm256_fmadd_ps(a, b, c); }
    inline real_v remainder( const real_v& a, const real_v& b ){ return a - real_v(_mm256_round_ps(a/b, _MM_FROUND_TO_NEG_INF)) * b; }
    inline real_v atan2( const real_v& y, const real_v& x ){ 
      const auto* bx = x.to_ptr();
      const auto* by = y.to_ptr();
      real_v rt; 
      for( unsigned i = 0 ; i != real_v::size ; ++i ) rt.to_ptr()[i] = std::atan2( by[i] , bx[i] );
      return rt; 
    }
    inline real_v gather( const double* base_addr, const real_v& offsets)
    {
      /// temporary -> improve by loading into two 256b SIMD registers, casting and merging 
      std::array<float, real_v::size> tmp;
      auto ptr = store( offsets.to_int() );
   //   int32_t* ptr = (int32_t*)(&ints); 
      for( int i = 0 ; i != real_v::size; ++i ) tmp[i] = real_v::scalar_type( base_addr[ptr[i]] );
      return real_v( tmp.data() ); 
    }
    
    inline real_v fmod( const real_v& a, const real_v& b )
    {
      auto r = remainder( abs(a), abs(b) );
      return select( a > 0., r, -r );
    }

    inline std::ostream& operator<<( std::ostream& os, const real_v& obj ) { 
      auto buffer = obj.to_array();
      for( unsigned i = 0 ; i != real_v::size; ++i ) os << buffer[i] << " ";
      return os; 
    }
 
    using complex_v = Complex<real_v>; 
    inline complex_v select(const real_v& mask, const complex_v& a, const complex_v& b ) { return complex_v( select(mask, a.real(), b.real()), select(mask, a.imag(), b.imag() ) ) ; }
    inline complex_v select(const real_v& mask, const real_v&   a, const complex_v& b ) { return complex_v( select(mask, a   , b.real()), select(mask, 0.f, b.imag()) ); }
    inline complex_v select(const real_v& mask, const complex_v& a, const real_v& b   ) { return complex_v( select(mask, a.real(), b )  , select(mask, a.imag(), 0.f) ); }
    inline complex_v select(const bool& mask   , const complex_v& a, const complex_v& b ) { return mask ? a : b; }
    #pragma omp declare reduction(+: real_v: \
	omp_out = omp_out + omp_in)
    #pragma omp declare reduction(+: complex_v: \
	omp_out = omp_out + omp_in)
  }  
}

#endif
