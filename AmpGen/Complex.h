#ifndef AMPGEN_COMPLEX_H 
#define AMPGEN_COMPLEX_H 1 

#include <complex>

namespace AmpGen {
  template <typename real_t>
  struct Complex {
    real_t re;
    real_t im;
    Complex() = default;
    Complex( const real_t& re, const real_t& im) : re(re), im(im) {}
    Complex( const float&   re, const float& im) : re(re), im(im) {}
    Complex( const std::complex<double>& f ) : re( f.real() ), im( f.imag() ) {}
    Complex( const std::complex<float>& f  ) : re( f.real() ), im( f.imag() ) {}
    explicit Complex( const real_t& arg ) : re(arg) {};
    explicit Complex( const double& arg ) : re(arg) {};
    inline Complex operator+=(const Complex& rhs ); 
    inline Complex operator-=(const Complex& rhs ); 
    inline Complex operator*=(const Complex& rhs ); 
    inline Complex operator/=(const Complex& rhs ); 
    real_t real() const { return re; }
    real_t imag() const { return im; }
    real_t norm() const { return re*re + im *im ; }
  };
  
  template<typename real_t> inline real_t            real(const Complex<real_t>& arg ){ return arg.re ; }
  template<typename real_t> inline real_t            imag(const Complex<real_t>& arg ){ return arg.im ; }
  template<typename real_t> inline real_t             abs(const Complex<real_t>& v ) { return sqrt( v.re * v.re + v.im * v.im ) ; }
  template<typename real_t> inline real_t            norm(const Complex<real_t>& v ) { return  ( v.re * v.re + v.im * v.im ) ; }
  template<typename real_t> inline Complex<real_t>   conj(const Complex<real_t>& arg ){ return Complex<real_t>(arg.re, -arg.im) ; }
  

  template<typename real_t> inline Complex<real_t> operator+( const Complex<real_t>& lhs, const real_t& rhs ) { return Complex<real_t>(lhs.re + rhs, lhs.im); }
  template<typename real_t> inline Complex<real_t> operator-( const Complex<real_t>& lhs, const real_t& rhs ) { return Complex<real_t>(lhs.re - rhs, lhs.im); }
  template<typename real_t> inline Complex<real_t> operator*( const Complex<real_t>& lhs, const real_t& rhs ) { return Complex<real_t>(lhs.re*rhs, lhs.im*rhs); }
  template<typename real_t> inline Complex<real_t> operator/( const Complex<real_t>& lhs, const real_t& rhs ) { return Complex<real_t>(lhs.re/rhs, lhs.im/rhs); }
  template<typename real_t> inline Complex<real_t> operator+( const real_t& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs + rhs.re,  rhs.im); }
  template<typename real_t> inline Complex<real_t> operator-( const real_t& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs - rhs.re, - rhs.im); }
  template<typename real_t> inline Complex<real_t> operator*( const real_t& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs*rhs.re, lhs*rhs.im); }
  template<typename real_t> inline Complex<real_t> operator/( const real_t& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>( lhs * rhs.re , -lhs *rhs.im) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
  template<typename real_t> inline Complex<real_t> operator+( const Complex<real_t>& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs.re + rhs.re, lhs.im + rhs.im); }
  template<typename real_t> inline Complex<real_t> operator-( const Complex<real_t>& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs.re - rhs.re, lhs.im - rhs.im); }
  template<typename real_t> inline Complex<real_t> operator*( const Complex<real_t>& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs.re*rhs.re - lhs.im*rhs.im, lhs.re*rhs.im  + lhs.im*rhs.re); }
  template<typename real_t> inline Complex<real_t> operator/( const Complex<real_t>& lhs, const Complex<real_t>& rhs ) { return Complex<real_t>(lhs.re*rhs.re + lhs.im*rhs.im, -lhs.re*rhs.im  + lhs.im*rhs.re) / (rhs.re * rhs.re + rhs.im * rhs.im ); }
  template<typename real_t> inline Complex<real_t> operator-( const Complex<real_t>& x ) { return -1.f * x; }
  template<typename real_t> inline Complex<real_t> Complex<real_t>::operator+=(const Complex<real_t>& rhs ){ *this = *this + rhs; return *this; }
  template<typename real_t> inline Complex<real_t> Complex<real_t>::operator-=(const Complex<real_t>& rhs ){ *this = *this - rhs; return *this; }
  template<typename real_t> inline Complex<real_t> Complex<real_t>::operator*=(const Complex<real_t>& rhs ){ *this = *this * rhs; return *this; }
  template<typename real_t> inline Complex<real_t> Complex<real_t>::operator/=(const Complex<real_t>& rhs ){ *this = *this / rhs; return *this; }
  template<typename real_t> inline Complex<real_t> exp( const Complex<real_t>& v ){
    auto [s,c] = sincos(v.im);
    return exp(v.re) * Complex<real_t>(c, s);
  }
  template <typename real_t> 
  inline Complex<real_t> sqrt( const Complex<real_t>& v )
  {
    auto r = abs(v);
    return Complex ( sqrt( 0.5 * (r + v.re) ), sign(v.im) * sqrt( 0.5*( r - v.re ) ) );
  }
  template <typename real_t>
  inline Complex<real_t> log( const Complex<real_t>& v )
  {
    return Complex<real_t>( 0.5 * log( v.norm() ) , atan2(v.im, v.re) );
  }
  template <typename real_t> 
  inline std::ostream& operator<<( std::ostream& os, const Complex<real_t>& obj ) { return os << "( "<< obj.re << ") (" << obj.im << ")"; }
}

#endif
