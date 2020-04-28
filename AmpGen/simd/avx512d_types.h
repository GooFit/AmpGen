#ifndef AMPGEN_AVXd_TYPES
#define AMPGEN_AVXd_TYPES 1

#include <immintrin.h>
#include <array>
#include <iostream>
#include <complex>
#include <omp.h>

namespace AmpGen {
  namespace AVX512d {
    #define stl_fallback( x ) \
      inline real_v x( const real_v& v ){  auto a = v.to_array(); return real_v( std::x(a[0]), std::x(a[1]), std::x(a[2]), std::x(a[3]), std::x(a[4]), std::x(a[5]), std::x(a[6]), std::x(a[7]) ) ; } 
       
    struct real_v {
      __m512d data;
      static constexpr unsigned size = 8;
      typedef double scalar_type;
      real_v() = default; 
      real_v(__m512d data ) : data(data) {}
      real_v(const double& f ) : data( _mm512_set1_pd( f )) {}
      real_v(
          const double& x0, const double& x1, const double& x2, const double& x3,
          const double& x4, const double& x5, const double& x6, const double& x7)
      {
        double tmp[8] = {x0,x1,x2,x3,x4,x5,x6,x7};
        data = _mm512_loadu_pd(tmp); 
      }
      real_v(const double* f ) : data( _mm512_loadu_pd( f ) ) {}
      void store( double* ptr ) const { _mm512_storeu_pd( ptr, data ); }
      std::array<double, 8> to_array() const { std::array<double, 8> b; store( &b[0] ); return b; }
      double at(const unsigned i) const { return to_array()[i] ; }       
      operator __m512d() const { return data ; } 
    };

    inline real_v operator+( const real_v& lhs, const real_v& rhs ) { return _mm512_add_pd(lhs, rhs); }
    inline real_v operator-( const real_v& lhs, const real_v& rhs ) { return _mm512_sub_pd(lhs, rhs); }
    inline real_v operator*( const real_v& lhs, const real_v& rhs ) { return _mm512_mul_pd(lhs, rhs); }
    inline real_v operator/( const real_v& lhs, const real_v& rhs ) { return _mm512_div_pd(lhs, rhs); }
    inline real_v operator-( const real_v& x ) { return -1.f * x; }
    inline real_v operator&( const real_v& lhs, const real_v& rhs ) { return _mm512_and_pd( lhs, rhs ); }
    inline real_v operator|( const real_v& lhs, const real_v& rhs ) { return _mm512_or_pd( lhs, rhs ); }
    inline real_v operator^( const real_v& lhs, const real_v& rhs ) { return _mm512_xor_pd( lhs, rhs ); }
    inline real_v operator+=(real_v& lhs, const real_v& rhs ){ lhs = lhs + rhs; return lhs; }
    inline real_v operator-=(real_v& lhs, const real_v& rhs ){ lhs = lhs - rhs; return lhs; }
    inline real_v operator*=(real_v& lhs, const real_v& rhs ){ lhs = lhs * rhs; return lhs; }
    inline real_v operator/=(real_v& lhs, const real_v& rhs ){ lhs = lhs / rhs; return lhs; }
    inline real_v operator&&( const real_v& lhs, const real_v& rhs ) { return _mm512_and_pd( lhs, rhs ); }
    inline real_v operator||( const real_v& lhs, const real_v& rhs ) { return _mm512_or_pd( lhs, rhs ); }
    inline real_v operator!( const real_v& x ) { return x ^ _mm512_castsi512_pd( _mm512_set1_epi32( -1 ) ); }
    inline __mmask8 operator<( const real_v& lhs, const real_v& rhs ) { return _mm512_cmp_pd_mask( lhs, rhs, _CMP_LT_OS ); }
    inline __mmask8 operator>( const real_v& lhs, const real_v& rhs ) { return _mm512_cmp_pd_mask( lhs, rhs, _CMP_GT_OS ); }
    inline __mmask8 operator==( const real_v& lhs, const real_v& rhs ){ return _mm512_cmp_pd_mask( lhs, rhs, _CMP_EQ_OS ); }
    inline real_v sqrt( const real_v& v ) { return _mm512_sqrt_pd(v); } 
    inline real_v abs ( const real_v& v ) { return _mm512_andnot_pd(_mm512_set1_pd(-0.), v);  }
    // inline real_v sin( const real_v& v )  { return sin512_pd(v) ; }
    // inline real_v cos( const real_v& v )  { return cos512_pd(v) ; }
    // inline real_v tan( const real_v& v )  { real_v s; real_v c; sincos512_pd(v, (__m512*)&s, (__m512*)&c) ; return s/c; }
    // inline real_v exp( const real_v& v )  { return exp512_ps(v) ; }
    inline real_v select(const __mmask8& mask, const real_v& a, const real_v& b ) { return _mm512_mask_mov_pd( b, mask, a ); }
    inline real_v select(const bool& mask   , const real_v& a, const real_v& b ) { return mask ? a : b; } 
    inline real_v sign  ( const real_v& v){ return select( v > 0., +1., -1. ); }
    inline real_v atan2( const real_v& y, const real_v& x ){ 
      std::array<double, 8> bx{x.to_array()}, by{y.to_array()};
      return real_v ( 
          std::atan2(by[0], bx[0]) , std::atan2( by[1], bx[1]), std::atan2( by[2], bx[2]), std::atan2( by[3], bx[3]) ,
          std::atan2(by[4], bx[4]) , std::atan2( by[5], bx[5]), std::atan2( by[6], bx[6]), std::atan2( by[7], bx[7]) );
    }
    inline __m512i double_to_int( const real_v& x )
    {
      auto xr = _mm512_roundscale_pd(x, _MM_FROUND_TO_ZERO);
      // based on:  https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx
      return _mm512_sub_epi64(_mm512_castpd_si512(_mm512_add_pd(xr, _mm512_set1_pd(0x0018000000000000))), 
             _mm512_castpd_si512(_mm512_set1_pd(0x0018000000000000)));
    }
    inline real_v gather( const double* base_addr, const real_v& offsets) 
    {
     return _mm512_i64gather_pd(double_to_int(offsets), base_addr, sizeof(double));
    } 
    
    inline void frexp(const real_v& value, real_v& mant, real_v& exponent)
    {
      auto arg_as_int = _mm512_castpd_si512(value);
      static const real_v offset(4503599627370496.0 + 1022.0);   // 2^52 + 1022.0
      static const __m512i pow2_52_i = _mm512_set1_epi64(0x4330000000000000); // *reinterpret_cast<const uint64_t*>(&pow2_52_d);
      auto b = _mm512_srl_epi64(arg_as_int, _mm_cvtsi32_si128(52));
      auto c = _mm512_or_si512( b , pow2_52_i);
      exponent = real_v( _mm512_castsi512_pd(c) ) - offset;
      mant     = _mm512_castsi512_pd(_mm512_or_si512(_mm512_and_si512 (arg_as_int, _mm512_set1_epi64(0x000FFFFFFFFFFFFFll) ), _mm512_set1_epi64(0x3FE0000000000000ll)));
    }
    
    inline real_v fmadd( const real_v& a, const real_v& b, const real_v& c )
    {
      return _mm512_fmadd_pd(a, b, c);
    }
    inline real_v log(const real_v& arg)
    {
      static const real_v corr     = 0.693147180559945286226764;
      static const real_v CL15     = 0.148197055177935105296783;
      static const real_v CL13     = 0.153108178020442575739679;
      static const real_v CL11     = 0.181837339521549679055568;
      static const real_v CL9      = 0.22222194152736701733275;
      static const real_v CL7      = 0.285714288030134544449368;
      static const real_v CL5      = 0.399999999989941956712869;
      static const real_v CL3      = 0.666666666666685503450651;
      static const real_v CL1      = 2.0;
      real_v mant, exponent;
      frexp(arg, mant, exponent);
      auto x  = (mant - 1.) / (mant + 1.);
      auto x2 = x * x;
      auto p = fmadd(CL15, x2, CL13);
      p = fmadd(p, x2, CL11);
      p = fmadd(p, x2, CL9);
      p = fmadd(p, x2, CL7);
      p = fmadd(p, x2, CL5);
      p = fmadd(p, x2, CL3);
      p = fmadd(p, x2, CL1);
      p = fmadd(p, x, corr * exponent);
      return p;
    }
    stl_fallback( exp )
    stl_fallback( tan )
    stl_fallback( sin )
    stl_fallback( cos )
    inline real_v remainder( const real_v& a, const real_v& b ){ return a - real_v(_mm512_roundscale_pd(a/b, _MM_FROUND_TO_NEG_INF)) * b; }
    inline real_v fmod( const real_v& a, const real_v& b )
    {
      auto r = remainder( abs(a), abs(b) );
      return select( a > 0., r, -r );
    }

    struct complex_v {
      real_v re;
      real_v im;
      static constexpr unsigned size = 8;
      typedef std::complex<double> scalar_type;

      real_v real() const { return re; }
      real_v imag() const { return im; }
      real_v norm() const { return re*re + im *im ; }
      complex_v() = default; 
      complex_v( const real_v& re, const real_v& im) : re(re), im(im) {}
      complex_v( const float&   re, const float& im) : re(re), im(im) {}
      complex_v( const std::complex<double>& f ) : re( f.real() ), im( f.imag() ) {}
      complex_v( const std::complex<float>& f  ) : re( f.real() ), im( f.imag() ) {}
      explicit complex_v( const real_v& arg ) : re(arg) {};
      explicit complex_v( const double& arg ) : re(arg) {};
      const std::complex<double> at(const unsigned i) const { return std::complex<float>(re.to_array()[i], im.to_array()[i]) ; }       
      void store( double* sre, double* sim ){ re.store(sre); im.store(sim);  } 
      void store( scalar_type* r ) const { 
        auto re_arr = re.to_array();
        auto im_arr = im.to_array();
        for( unsigned i = 0 ; i != re_arr.size(); ++i ) r[i] = scalar_type( re_arr[i], im_arr[i] ); 
      }
      auto to_array() const 
      {
        std::array<scalar_type, 8> rt;
        store( rt.data() );
        return rt;
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
    inline complex_v select(const __mmask8& mask, const complex_v& a, const complex_v& b ) { return complex_v( select(mask, a.re, b.re), select(mask, a.im, b.im ) ) ; }
    inline complex_v select(const __mmask8& mask, const real_v&   a, const complex_v& b ) { return complex_v( select(mask, a   , b.re), select(mask, 0.f, b.im) ); }
    inline complex_v select(const __mmask8& mask, const complex_v& a, const real_v& b   ) { return complex_v( select(mask, a.re, b )  , select(mask, a.im, 0.f) ); }
    inline complex_v select(const bool& mask   , const complex_v& a, const complex_v& b ) { return mask ? a : b; }
    inline complex_v exp( const complex_v& v ){ 
      return exp( v.re) * complex_v( cos( v.im ), sin( v.im ) );
    }
    inline complex_v sqrt( const complex_v& v )
    {
      auto r = abs(v);
      return complex_v ( sqrt( 0.5 * (r + v.re) ), sign(v.im) * sqrt( 0.5*( r - v.re ) ) );
    }
    inline complex_v log( const complex_v& v )
    {
      return complex_v( log( v.re ) , atan2(v.im, v.re) );
    } 

    inline std::ostream& operator<<( std::ostream& os, const complex_v& obj ) { return os << "( "<< obj.re << ") (" << obj.im << ")"; }
    #pragma omp declare reduction(+: real_v: \
	omp_out = omp_out + omp_in) 
    #pragma omp declare reduction(+: complex_v: \
	omp_out = omp_out + omp_in) 
  
  }  
}

#endif
