#include "AmpGen/Wigner.h"
#include "AmpGen/Simplify.h"

using namespace AmpGen;

double factorial( const double& z )
{
  double f=1;
  for( int i=1;i<=z;++i) f*=i;
  return f;
}

double nCr( const int& n, const int& r ){
  double z=1;
  for( int f=1; f <= r ; ++f ) z *= double(n+1-f)/double(f);
  return z;
}

/// computes (1+x)^n
Expression ExpandedBinomial( const Expression& x, const unsigned int& n ){
  Expression sum;
  for( unsigned int k = 0 ; k <= n ; ++k ) sum = sum + nCr(n,k) * fcn::fpow(x,k);
  return sum; 
}


Expression AmpGen::wigner_d( const Expression& beta, const double& j, const double& m, const double& n )
{
  int k_min = std::max(0.,m+n);
  int k_max = std::min(j+m,j+n);
  Expression cb = fcn::cos(beta);
  Expression sgn_sin = Ternary(beta>0, 1, -1 );
  Expression sum = 0 ; 
  double w2_num = factorial(j+m) * factorial(j-m) * factorial(j+n) * factorial(j-n);
  double nc_intpart = 0;
  double ns_intpart = 0; 
  double frac_nc = modf( k_min -(m+n)/2.    , &nc_intpart );
  double frac_ns = modf( j + (m+n)/2. -k_min, &ns_intpart );
  Expression fractional_part = 1 ; 
  Expression sign = pow(-1,j+m) * fcn::pow( sgn_sin, frac_ns==0.5);
  
  if( frac_nc == 0.5 && frac_ns != 0.5 )      fractional_part = fcn::sqrt(1+cb);
  else if( frac_nc != 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb);
  else if( frac_nc == 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb*cb);  
  for( double k = k_min; k <= k_max ; ++k )
  {
    double w_den  = factorial(k) * factorial(j+m-k)*factorial(j+n-k) * factorial(k-m-n);
    double norm = pow(-1,k) * sqrt(w2_num)/( w_den * pow(2,j) );
    Expression p1 = ExpandedBinomial( cb, int( k - (m+n)/2.)    );
    Expression p2 = ExpandedBinomial(-cb, int( j + (m+n)/2. - k)); 
    sum = sum + norm * p1 * p2;  
  }
  auto simplified = NormalOrderedExpression(sum);
  return sign * fractional_part * simplified; 
}
