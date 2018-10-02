#include "AmpGen/Wigner.h"
#include "AmpGen/Simplify.h"

using namespace AmpGen;

double fact( const double& z )
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


Expression AmpGen::wigner_d( const Expression& cb, const double& j, const double& m, const double& n )
{
  int k_min = std::max(0.,m+n);
  int k_max = std::min(j+m,j+n);
  //Expression cb = fcn::cos(beta);
  //Expression sgn_sin = Ternary(beta>0, 1, -1 );
  Expression sum = 0 ; 
  double w2_num = fact(j+m) * fact(j-m) * fact(j+n) * fact(j-n);
  double nc_intpart = 0;
  double ns_intpart = 0; 
  double frac_nc = modf( k_min -(m+n)/2.    , &nc_intpart );
  double frac_ns = modf( j + (m+n)/2. -k_min, &ns_intpart );
  Expression fractional_part = 1 ; 
  //Expression sign = pow(-1,j+m) * fcn::pow( sgn_sin, frac_ns==0.5);

  if( frac_nc == 0.5 && frac_ns != 0.5 )      fractional_part = fcn::sqrt(1+cb);
  else if( frac_nc != 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb);
  else if( frac_nc == 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb*cb);  
  for( double k = k_min; k <= k_max ; ++k )
  {
    double w_den  = fact(k) * fact(j+m-k)*fact(j+n-k) * fact(k-m-n);
    double norm = pow(-1,k) * sqrt(w2_num)/( w_den * pow(2,j) );
    Expression p1 = ExpandedBinomial( cb, int( k - (m+n)/2.)    );
    Expression p2 = ExpandedBinomial(-cb, int( j + (m+n)/2. - k)); 
    sum = sum + norm * p1 * p2;  
  }
  auto simplified = NormalOrderedExpression(sum);
  return fractional_part * simplified; 
}

double AmpGen::CG( 
    const double& j1,
    const double& m1,
    const double& j2, 
    const double& m2,
    const double& J,
    const double& M )
{
  if( m1+m2!=M ) return 0;
  double f1 = (2*J+1)*fact(J+j1-j2)*fact(J-j1+j2)*fact(j1+j2-J) ;
  double f2 = fact(j1+m1)*fact(j1-m1)*fact(j2+m2)*fact(j2-m2)*fact(J+M)*fact(J-M);

  double norm = f1 * f2 / fact(J+j1+j2+1) ;
  double sum  = 0;

  for( int nu=0; nu <= j1+j2-J ; ++nu){
    double arg1 = j1+j2-J-double(nu);
    double arg2 = j1  -m1-double(nu);
    double arg3 = j2+  m2-double(nu);
    double arg4 = J-j2+m1+double(nu);
    double arg5 = J-j1-m2+double(nu);
    if( arg1 < 0 || arg2 < 0 || arg3 < 0 || arg4 < 0 || arg5 < 0 ) continue ; 
    int sgn = nu % 2 == 0 ? 1 : -1;
    double to_add =  sgn / (fact(nu)*fact(arg1)*fact(arg2)*fact(arg3)*fact(arg4)*fact(arg5) );
    sum = sum + to_add ; 
  }
  return sqrt(norm) * sum ; 
}

using namespace AmpGen::fcn;

Tensor AmpGen::rotationMatrix( const Tensor& P ){
  Tensor R(std::vector<size_t>({4,4}));
  auto px = P[0];
  auto py = P[1];
  auto pz = P[2];
  auto p2  = make_cse(px*px + py*py + pz*pz);
  auto pt2 = make_cse(px*px + py*py);
  auto p   = sqrt(p2); 
  Expression f    = Ternary( abs(p  ) < 1e-8  , 0 , pz/p - 1);
  Expression pxn  = Ternary( abs(pt2) < 1e-8  , 0 , px/sqrt(pt2));
  Expression pyn  = Ternary( abs(pt2) < 1e-8  , 0 , py/sqrt(pt2));
  Expression pxn2 = Ternary( abs(px)  < 1e-8  , 0 , px/p);
  Expression pyn2 = Ternary( abs(py)  < 1e-8  , 0 , py/p);
  R(0,0) = 1. + pxn*pxn*f;
  R(1,0) = pxn*pyn*f;
  R(2,0) = pxn2;
  R(0,1) = pxn*pyn*f;
  R(1,1) = 1. + pyn*pyn*f;
  R(2,1) = pyn2;
  R(0,2) = -pxn2;
  R(1,2) = -pyn2;
  R(2,2) = 1+f; 
  R(3,3) = 1.0;
  return R;
}

Tensor AmpGen::helicityTransformMatrix( const Tensor& P, const Expression& E, const Expression& M, const double& ve )
{
  auto dim = std::vector<size_t>({4,4});
  Tensor L(dim);
  Tensor R = rotationMatrix(P);
  Tensor::Index a,b,c;
  if( ve == -1 ){
    Tensor Rx(dim);
    Rx(0,0) =  1.0;
    Rx(1,1) = -1.0;
    Rx(2,2) = -1.0;
    Rx(3,3) =  1.0;
    R = Rx(a,b) * R(b,c);
  }
  Expression p = sqrt( make_cse( P[0]*P[0] + P[1]*P[1] + P[2]*P[2] ) );
  L(0,0) = 1.;
  L(1,1) = 1.;
  L(2,2) =  E / M;
  L(2,3) =  - p / M;
  L(3,2) =  - p / M;
  L(3,3) =  E / M;
  return L(a,b) * R(b,c);
}

Expression AmpGen::wigner_D(const Tensor& P, 
                            const double& J, 
                            const double& lA, 
                            const double& lB, 
                            const double& lC )
{
  Expression pz = make_cse( P[2] / sqrt( P[0]*P[0] + P[1] * P[1] + P[2]*P[2] ) );  
  Expression pt2 = make_cse( P[0]*P[0] + P[1]*P[1] );
  Expression px = P[0] / sqrt( pt2 );
  Expression py = P[1] / sqrt( pt2 );
  Expression I(std::complex<double>(0,1));
  auto little_d = make_cse ( wigner_d( pz, J, lA, lB-lC ) );
  return  fpow( px - I * py, lB-lC-lA) * little_d; 
}

