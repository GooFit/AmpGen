#ifndef AMPGEN_AVX_TYPES
#define AMPGEN_AVX_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#include <omp.h>

#if USE_MVEC
extern "C" void    _ZGVdN8vvv_sincos(__m256 x, __m256i ptrs, __m256i ptrc);
#endif

#if USE_MVEC 
#define libmvec_alias( function_name) \
  extern "C" __m256 _ZGVcN8v_##function_name(__m256 x);                                    \
  inline real_v function_name( const real_v& v ){ return _ZGVcN8v_##function_name (v) ; }
#else
#define libmvec_alias( F ) \
  inline real_v F( const real_v& v ){ auto arr = v.to_ptr(); return real_v(  \
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
      real_v(const float& f )      : data( _mm256_set1_ps(f) ) {}
      real_v(const double& f )     : data( _mm256_set1_ps( float(f) )) {}
      explicit real_v(const float* f )      : data( _mm256_loadu_ps( f ) ) {}
      real_v(const float& x0, const float& x1, const float& x2, const float& x3,
          const float& x4, const float& x5, const float& x6, const float& x7)
      {
        data = _mm256_set_ps(x7,x6,x5,x4,x3,x2,x1,x0); 
      }

      void store( float* ptr ) const { _mm256_storeu_ps( ptr, data ); }
      std::array<float, 8> to_array() const { std::array<float, 8> b; store( &b[0] ); return b; }
      const float* to_ptr() const { return reinterpret_cast<const double*>( &data ) ; }
      float at(const unsigned i) const { return to_ptr()[i] ; }       
      operator __m256() const { return data ; } 
      inline real_v operator+=(const real_v& rhs ){ *this = *this + rhs; return *this; }
      inline real_v operator-=(const real_v& rhs ){ *this = *this - rhs; return *this; }
      inline real_v operator*=(const real_v& rhs ){ *this = *this * rhs; return *this; }
      inline real_v operator/=(const real_v& rhs ){ *this = *this / rhs; return *this; }
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
    inline real_v operator<( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_LE_OS ); }
    inline real_v operator>( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_ps( lhs, rhs, _CMP_GE_OS ); }
    inline real_v operator==( const real_v& lhs, const real_v& rhs ){ return _mm256_cmp_ps( lhs, rhs, _CMP_EQ_OS ); }
    inline real_v sqrt( const real_v& v ) { return _mm256_sqrt_ps(v); } 
    libmvec_alias(sin)
    libmvec_alias(cos)
    libmvec_alias(exp)
    libmvec_alias(log)
    
    inline void sincos( const real_v& v, real_v& s, real_v& c )
    {
#if USE_MVEC
      __m256i sp = _mm256_add_epi64(_mm256_set1_epi32x((uint64_t)&s),_mm256_set_epi32x(28,24,20,16,12,8,4,0)); 
      __m256i cp = _mm256_add_epi64(_mm256_set1_epi32x((uint64_t)&c),_mm256_set_epi32x(28,24,20,16,12,8,4,0)); 
      _ZGVdN8vvv_sincos(v,sp,cp); 
#else
      s = sin(v);
      c = cos(v);
#endif
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
    
    inline real_v abs ( const real_v& v ) { return v & _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ); }
    inline real_v select(const real_v& mask, const real_v& a, const real_v& b ) { return _mm256_blendv_ps( b, a, mask ); }
    inline real_v select(const bool& mask   , const real_v& a, const real_v& b ) { return mask ? a : b; } 
    inline real_v atan2( const real_v& y, const real_v& x ){ 
      const auto* bx = reinterpret_cast<const double*>( &x.data );
      const auto* by = reinterpret_cast<const double*>( &y.data ); 
      for( unsigned i = 0 ; i != 8 ; ++i ) rt[i] = std::atan2( by[i] , bx[i] );
      return real_v (rt.data() ); 
    }
    
    inline real_v fmadd( const real_v& a, const real_v& b, const real_v& c )
    {
      return _mm256_fmadd_ps(a, b, c);
    }
    inline real_v remainder( const real_v& a, const real_v& b ){ return a - real_v(_mm256_round_ps(a/b, _MM_FROUND_TO_NEG_INF)) * b; }
    inline real_v fmod( const real_v& a, const real_v& b )
    {
      auto r = remainder( abs(a), abs(b) );
      return select( a > 0., r, -r );
    }

    inline std::ostream& operator<<( std::ostream& os, const real_v& obj ) { 
      auto buffer = obj.to_array();
      for( unsigned i = 0 ; i != 8; ++i ) os << buffer[i] << " ";
      return os; 
    }
 
    using complex_v = std::complex<real_v>; 
    inline complex_v operator+( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.real() + rhs, lhs.imag()); }
    inline complex_v operator-( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.real() - rhs, lhs.imag()); }
    inline complex_v operator*( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.real()*rhs, lhs.imag()*rhs); }
    inline complex_v operator/( const complex_v& lhs, const real_v& rhs ) { return complex_v(lhs.real()/rhs, lhs.imag()/rhs); }
    inline complex_v operator+( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs + rhs.real(),  rhs.imag()); }
    inline complex_v operator-( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs - rhs.real(), - rhs.imag()); }
    inline complex_v operator*( const real_v& lhs, const complex_v& rhs ) { return complex_v(lhs*rhs.real(), lhs*rhs.imag()); }
    inline complex_v operator/( const real_v& lhs, const complex_v& rhs ) { return complex_v( lhs * rhs.real() , -lhs *rhs.imag()) / (rhs.real() * rhs.real() + rhs.imag() * rhs.imag() ); }
    inline real_v abs( const complex_v& v ) { return sqrt( v.real() * v.real() + v.imag() * v.imag() ) ; }
    inline real_v norm( const complex_v& v ) { return  ( v.real() * v.real() + v.imag() * v.imag() ) ; }
    inline complex_v select(const real_v& mask, const complex_v& a, const complex_v& b ) { return complex_v( select(mask, a.real(), b.real()), select(mask, a.imag(), b.imag() ) ) ; }
    inline complex_v select(const real_v& mask, const real_v&   a, const complex_v& b ) { return complex_v( select(mask, a   , b.real()), select(mask, 0.f, b.imag()) ); }
    inline complex_v select(const real_v& mask, const complex_v& a, const real_v& b   ) { return complex_v( select(mask, a.real(), b )  , select(mask, a.imag(), 0.f) ); }
    inline complex_v select(const bool& mask   , const complex_v& a, const complex_v& b ) { return mask ? a : b; }
    inline complex_v exp( const complex_v& v ){
      auto [s,c] = sincos( v.imag());
      return exp(v.real()) * complex_v(c, s);
    }
    inline complex_v sqrt( const complex_v& v )
    {
      auto r = abs(v);
      return complex_v ( sqrt( 0.5 * (r + v.real()) ), sign(v.imag()) * sqrt( 0.5*( r - v.real() ) ) );
    }
    inline complex_v log( const complex_v& v )
    {
      return complex_v( 0.5 * log( norm(v) ) , atan2(v.imag(), v.real()) );
    }

    inline std::ostream& operator<<( std::ostream& os, const complex_v& obj ) { return os << "( "<< obj.real() << ") (" << obj.imag() << ")"; }
    #pragma omp declare reduction(+: real_v: \
	omp_out = omp_out + omp_in)
    #pragma omp declare reduction(+: complex_v: \
	omp_out = omp_out + omp_in)
  }  
}

#endif
