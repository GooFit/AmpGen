#ifndef AMPGEN_ARM128d_TYPES
#define AMPGEN_ARM128d_TYPES 1

#include <arm_neon.h>
#include <array>
#include <iostream>
#include <complex>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>


#if USE_MVEC
extern "C" void    _ZGVdN4vvv_sincos(__m256d x, __m256i ptrs, __m256i ptrc);
#define libmvec_alias( function_name) \
  extern "C" __m256d _ZGVcN4v_##function_name(__m256d x);                                    \
inline real_v function_name( const real_v& v ){ return _ZGVcN4v_##function_name (v) ; }
#else
#define libmvec_alias( F ) \
  inline real_v F( const real_v& v ){ auto arr = v.to_ptr(); return real_v( std::F(arr[0]), std::F(arr[1]) ); }
#endif

namespace AmpGen {
  namespace ARM128d {

    struct real_v {
      float64x2_t data;
      static constexpr unsigned size = 2;
      typedef double scalar_type;
      real_v() = default;
      real_v(float64x2_t data ) : data(data) {}
      explicit real_v( uint64x2_t&& data ) : data( vcvtq_f64_u64(data) ) {}
      explicit real_v(  int64x2_t&& data ) : data( vcvtq_f64_s64(data) ) {}
      real_v(const scalar_type& f ) : data(vmovq_n_f64(f)) {}
      real_v(const scalar_type& x0, const scalar_type& x1 )
      {
        data = vsetq_lane_f64(x0, data, 0 );
        data = vsetq_lane_f64(x1, data, 1 ); 
      }
      explicit real_v(const scalar_type* f ) : data(  vld1q_f64( f ) ) {}
      real_v(const std::array<scalar_type, size> f ) : data(  vld1q_f64( f.data() ) ) {}
      void store( scalar_type* ptr ) const { vst1q_f64(ptr, data); }
      const scalar_type* to_ptr() const { return reinterpret_cast<const scalar_type*>( &data ) ; }
      scalar_type* to_ptr()       { return reinterpret_cast<scalar_type*>( &data ) ; }
      std::array<scalar_type, size> to_array() const { std::array<scalar_type, size> b; store( &b[0] ); return b; }
      int64x2_t to_int() const { return vcvtq_s64_f64(data); }
      double at(const unsigned i) const { return to_ptr()[i]; }
      operator float64x2_t() const { return data ; }
      inline real_v operator+=(const real_v& rhs );
      inline real_v operator-=(const real_v& rhs );
      inline real_v operator*=(const real_v& rhs );
      inline real_v operator/=(const real_v& rhs );
    };
    struct int_v { int_v( uint64x2_t&& data) : data(data){}; uint64x2_t data; operator uint64x2_t() const { return data;} };


    inline real_v operator+( const real_v& lhs, const real_v& rhs ) { return  vaddq_f64(lhs, rhs); }
    inline real_v operator-( const real_v& lhs, const real_v& rhs ) { return  vsubq_f64(lhs, rhs); }
    inline real_v operator*( const real_v& lhs, const real_v& rhs ) { return  vmulq_f64(lhs, rhs); }
    inline real_v operator/( const real_v& lhs, const real_v& rhs ) { return  vdivq_f64(lhs, rhs); }
    inline real_v operator-( const real_v& x ) { return -1.f * x; }

    // inline real_v operator&( const real_v& lhs, const real_v& rhs ) { return real_v( vceqq_f64( lhs, rhs ) ); }
    // inline real_v operator|( const real_v& lhs, const real_v& rhs ) { return _mm256_or_pd( lhs, rhs ); }
    // inline real_v operator^( const real_v& lhs, const real_v& rhs ) { return _mm256_xor_pd( lhs, rhs ); }
    inline int_v operator&&( const int_v& lhs, const int_v& rhs ) { return vandq_u64( lhs, rhs ); }
    inline int_v operator||( const int_v& lhs, const int_v& rhs ) { return vorrq_u64( lhs, rhs ); }
    // inline real_v operator!( const real_v& x ) { return x ^ _mm256_castsi256_pd( _mm256_set1_epi32( -1 ) ); }

    inline int_v operator<( const real_v& lhs, const real_v& rhs ) { return vcltq_f64(lhs,rhs); }
    inline int_v operator>( const real_v& lhs, const real_v& rhs ) { return vcgtq_f64(lhs,rhs); }
    inline int_v operator<=( const real_v& lhs, const real_v& rhs ){ return vcleq_f64( lhs, rhs ); }
    inline int_v operator>=( const real_v& lhs, const real_v& rhs ){ return vcleq_f64( lhs, rhs ); }
    inline int_v operator==( const real_v& lhs, const real_v& rhs ){ return vceqq_f64( lhs, rhs); }
    inline real_v sqrt( const real_v& v ) { return vsqrtq_f64(v); }
    inline real_v abs ( const real_v& v ) { return vabsq_f64(v); } 
    inline real_v real_v::operator+=(const real_v& rhs ){ *this = *this + rhs; return *this; }
    inline real_v real_v::operator-=(const real_v& rhs ){ *this = *this - rhs; return *this; }
    inline real_v real_v::operator*=(const real_v& rhs ){ *this = *this * rhs; return *this; }
    inline real_v real_v::operator/=(const real_v& rhs ){ *this = *this / rhs; return *this; }
    libmvec_alias( sin )
      libmvec_alias( cos )
      libmvec_alias( exp )
      libmvec_alias( log )
      inline void sincos( const real_v& v, real_v& s, real_v& c )
      {
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
    inline std::array<uint64_t, real_v::size> store( const int_v& v )
    {
      std::array<uint64_t, real_v::size> rt;
      vst1q_u64( rt.data(), v );
      return rt; 
    }

    inline real_v select(const int_v& mask, const real_v& a, const real_v& b ) { return vbslq_f64(mask, a, b); }
    inline real_v select(const bool& mask  , const real_v& a, const real_v& b ) { return mask ? a : b; }
    inline real_v sign  ( const real_v& v){ return select( v > real_v(0.), +1., -1. ); }
    inline real_v atan2( const real_v& y, const real_v& x ){
      const double* bx = x.to_ptr();
      const double* by = y.to_ptr(); 
      return real_v ( std::atan2(by[0], bx[0]), std::atan2( by[1], bx[1]) ); 
    }
    inline real_v gather( const double* base_addr, const real_v& offsets)
    {
      std::array<int64_t, real_v::size> offsets_p; 
      vst1q_s64( offsets_p.data(), offsets.to_int() ); 
      return real_v ( base_addr[offsets_p[0]], base_addr[offsets_p[1]] ); 
    }
    inline real_v fmadd( const real_v& a, const real_v& b, const real_v& c )
    {
      return vmlaq_f64(a, b, c);
    }
    inline real_v remainder( const real_v& a, const real_v& b ){ return a - b * real_v(vcvtq_u64_f64(a/b)); }
    inline real_v fmod( const real_v& a, const real_v& b ){ return remainder( abs(a), abs(b) ) * sign(a); }

    struct complex_v {
      real_v re;
      real_v im;
      typedef std::complex<double> scalar_type;
      static constexpr unsigned size = real_v::size;

      real_v real() const { return re; }
      real_v imag() const { return im; }
      real_v norm() const { return re*re + im *im ; }
      complex_v() = default;
      complex_v( const real_v& re, const real_v& im) : re(re), im(im) {}
      complex_v( const float&   re, const float& im) : re(re), im(im) {}
      complex_v( const std::complex<double>& f ) : re( f.real() ), im( f.imag() ) {}
      complex_v( const std::complex<float>& f  ) : re( f.real() ), im( f.imag() ) {}
      complex_v( const std::complex<double>* arr ) :
        re ( arr[0].real(), arr[1].real()),
        im ( arr[0].imag(), arr[1].imag()){}
      explicit complex_v( const real_v& arg ) : re(arg) {};
      explicit complex_v( const double& arg ) : re(arg) {};
      const std::complex<double> at(const unsigned i) const { return std::complex<double>(re.at(i), im.at(i)); }
      void store( double* sre, double* sim ){ re.store(sre); im.store(sim);  }
      void store( scalar_type* r ) const {
        auto re_arr = re.to_array();
        auto im_arr = im.to_array();
        for( unsigned i = 0 ; i != re_arr.size(); ++i ) r[i] = scalar_type( re_arr[i], im_arr[i] );
      }
      auto to_array() const
      {
        std::array<scalar_type, 2> rt;
        store( rt.data() );
        return rt;
      }
      inline complex_v operator+=(const complex_v& rhs ); 
      inline complex_v operator-=(const complex_v& rhs ); 
      inline complex_v operator*=(const complex_v& rhs ); 
      inline complex_v operator/=(const complex_v& rhs ); 
    };

    inline std::ostream& operator<<( std::ostream& os, const real_v& obj ) {
      auto data = obj.to_ptr();
      for( unsigned i = 0 ; i != 4; ++i ) os << data[i] << " ";
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
    inline complex_v select(const int_v& mask, const complex_v& a, const complex_v& b ) { return complex_v( select(mask, a.re, b.re), select(mask, a.im, b.im ) ) ; }
    inline complex_v select(const int_v& mask, const real_v&   a, const complex_v& b ) { return complex_v( select(mask, a   , b.re), select(mask, 0.f, b.im) ); }
    inline complex_v select(const int_v& mask, const complex_v& a, const real_v& b   ) { return complex_v( select(mask, a.re, b )  , select(mask, a.im, 0.f) ); }
    inline complex_v select(const bool& mask   , const complex_v& a, const complex_v& b ) { return mask ? a : b; }
    inline complex_v complex_v::operator+=(const complex_v& rhs ){ *this = *this + rhs; return *this; }
    inline complex_v complex_v::operator-=(const complex_v& rhs ){ *this = *this - rhs; return *this; }
    inline complex_v complex_v::operator*=(const complex_v& rhs ){ *this = *this * rhs; return *this; }
    inline complex_v complex_v::operator/=(const complex_v& rhs ){ *this = *this / rhs; return *this; }
    inline complex_v exp( const complex_v& v ){
      auto [s,c] = sincos( v.im);
      return exp(v.re) * complex_v(c, s);
    }
    inline complex_v sqrt( const complex_v& v )
    {
      auto r = abs(v);
      return complex_v ( sqrt( 0.5 * (r + v.re) ), sign(v.im) * sqrt( 0.5*( r - v.re ) ) );
    }
    inline complex_v log( const complex_v& v )
    {
      return complex_v( 0.5 * log( v.norm() ) , atan2(v.im, v.re) );
    }

    inline std::ostream& operator<<( std::ostream& os, const complex_v& obj ) { return os << "( "<< obj.re << ") (" << obj.im << ")"; }

  }
}

#endif
