#ifndef AMPGEN_AVXd_TYPES
#define AMPGEN_AVXd_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>

#if USE_MVEC
extern "C" void    _ZGVdN4vvv_sincos(__m256d x, __m256i ptrs, __m256i ptrc);
#endif

#if USE_MVEC
#define libmvec_alias( function_name) \
  extern "C" __m256d _ZGVcN4v_##function_name(__m256d x);                                    \
  inline real_v function_name( const real_v& v ){ return _ZGVcN4v_##function_name (v) ; }
#else
#define libmvec_alias( F ) \
  inline real_v F( const real_v& v ){ auto arr = v.to_ptr(); return real_v( std::F(arr[0]), std::F(arr[1]), std::F(arr[2]), std::F(arr[3])) ; }
#endif

namespace AmpGen {
  namespace AVX2d {

    struct real_v {
      __m256d data;
      static constexpr unsigned size = 4;
      typedef double scalar_type;
      real_v() = default;
      real_v(__m256d data ) : data(data) {}
      real_v(const double& f ) : data( _mm256_set1_pd( f )) {}
      real_v(const double& x0, const double& x1, const double& x2, const double& x3 )
      {
        data = _mm256_set_pd(x3,x2,x1,x0);
      }
      explicit real_v(const double* f ) : data( _mm256_loadu_pd( f ) ) {}
      real_v(const std::array<double,4> f ) : data( _mm256_loadu_pd( f.data() ) ) {}
      void store( double* ptr ) const { _mm256_storeu_pd( ptr, data ); }
      const double* to_ptr() const { return reinterpret_cast<const double*>( &data ) ; }
            double* to_ptr()       { return reinterpret_cast<double*>( &data ) ; }
      std::array<double, 4> to_array() const { std::array<double, 4> b; store( &b[0] ); return b; }
      double at(const unsigned i) const { return to_ptr()[i]; }
      operator __m256d() const { return data ; }
      inline real_v operator+=(const real_v& rhs );
      inline real_v operator-=(const real_v& rhs );
      inline real_v operator*=(const real_v& rhs );
      inline real_v operator/=(const real_v& rhs );
      inline __m256i to_int() const
      {
        // based on:  https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
        auto xr = _mm256_round_pd(data, _MM_FROUND_TO_NEG_INF);
        return _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(xr, _mm256_set1_pd(0x0018000000000000))),
               _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000)));
      }
    };

    inline real_v operator+( const real_v& lhs, const real_v& rhs ) { return _mm256_add_pd(lhs, rhs); }
    inline real_v operator-( const real_v& lhs, const real_v& rhs ) { return _mm256_sub_pd(lhs, rhs); }
    inline real_v operator*( const real_v& lhs, const real_v& rhs ) { return _mm256_mul_pd(lhs, rhs); }
    inline real_v operator/( const real_v& lhs, const real_v& rhs ) { return _mm256_div_pd(lhs, rhs); }
    inline real_v operator-( const real_v& x ) { return -1.f * x; }
    inline real_v operator&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_pd( lhs, rhs ); }
    inline real_v operator|( const real_v& lhs, const real_v& rhs ) { return _mm256_or_pd( lhs, rhs ); }
    inline real_v operator^( const real_v& lhs, const real_v& rhs ) { return _mm256_xor_pd( lhs, rhs ); }
    inline real_v operator&&( const real_v& lhs, const real_v& rhs ) { return _mm256_and_pd( lhs, rhs ); }
    inline real_v operator||( const real_v& lhs, const real_v& rhs ) { return _mm256_or_pd( lhs, rhs ); }
    inline real_v operator!( const real_v& x ) { return x ^ _mm256_castsi256_pd( _mm256_set1_epi32( -1 ) ); }
    inline real_v operator<( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_LT_OS ); }
    inline real_v operator>( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_GT_OS ); }
    inline real_v operator<=( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_LE_OS ); }
    inline real_v operator>=( const real_v& lhs, const real_v& rhs ) { return _mm256_cmp_pd( lhs, rhs, _CMP_GE_OS ); }
    inline real_v operator==( const real_v& lhs, const real_v& rhs ){ return _mm256_cmp_pd( lhs, rhs, _CMP_EQ_OS ); }
    inline real_v sqrt( const real_v& v ) { return _mm256_sqrt_pd(v); }
    inline real_v abs ( const real_v& v ) { return _mm256_andnot_pd(_mm256_set1_pd(-0.), v);  }
    inline real_v real_v::operator+=(const real_v& rhs ){ *this = *this + rhs; return *this; }
    inline real_v real_v::operator-=(const real_v& rhs ){ *this = *this - rhs; return *this; }
    inline real_v real_v::operator*=(const real_v& rhs ){ *this = *this * rhs; return *this; }
    inline real_v real_v::operator/=(const real_v& rhs ){ *this = *this / rhs; return *this; }
    libmvec_alias( sin )
    libmvec_alias( cos )
    libmvec_alias( exp )
    libmvec_alias( log )
    //inline real_v log( const real_v& v ){ auto arr = v.to_array(); return real_v( std::log(arr[0]), std::log(arr[1]), std::log(arr[2]), std::log(arr[3])) ; }
    inline void sincos( const real_v& v, real_v& s, real_v& c )
    {
#if USE_MVEC
      __m256i sp = _mm256_add_epi64(_mm256_set1_epi64x((uint64_t)&s),_mm256_set_epi64x(24,16,8,0));
      __m256i cp = _mm256_add_epi64(_mm256_set1_epi64x((uint64_t)&c),_mm256_set_epi64x(24,16,8,0));
      _ZGVdN4vvv_sincos(v,sp,cp);
#else
      s = sin(v);
      c = cos(v);
#endif
    }
    inline std::array<int64_t, real_v::size> store( const __m256i& v )
    {
      alignas(32) std::array<int64_t, real_v::size> rt; 
       _mm256_store_si256( (__m256i*)&rt[0], v);
       return rt; 
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

    inline real_v select(const real_v& mask, const real_v& a, const real_v& b ) { return _mm256_blendv_pd( b, a, mask ); }
    inline real_v select(const bool& mask   , const real_v& a, const real_v& b ) { return mask ? a : b; }
    inline real_v sign  ( const real_v& v){ return select( v > 0., +1., -1. ); }
    inline real_v conj( const real_v& v ) { return v; }
    
    inline real_v atan2( const real_v& y, const real_v& x ){
      const double* bx = x.to_ptr();
      const double* by = y.to_ptr(); 
      return real_v (std::atan2(by[0], bx[0]), std::atan2( by[1], bx[1]), std::atan2( by[2], bx[2]), std::atan2( by[3], bx[3]) );
    }
    inline real_v gather( const double* base_addr, const real_v& offsets)
    {
      return _mm256_i64gather_pd(base_addr, offsets.to_int(),sizeof(double));
    }
    inline real_v fmadd( const real_v& a, const real_v& b, const real_v& c )
    {
      return _mm256_fmadd_pd(a, b, c);
    }
    inline real_v remainder( const real_v& a, const real_v& b ){ return a - real_v(_mm256_round_pd(a/b, _MM_FROUND_TO_NEG_INF)) * b; }
    inline real_v fmod( const real_v& a, const real_v& b )
    {
      auto r = remainder( abs(a), abs(b) );
      return select( a > 0., r, -r );
    }

    inline std::ostream& operator<<( std::ostream& os, const real_v& obj ) {
      auto data = obj.to_ptr();
      for( unsigned i = 0 ; i != 4; ++i ) os << data[i] << " ";
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
