#include <nppdefs.h>
#include <iostream>
#include <math.h>
#include <array>
#include <cuComplex.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace ampgen_cuda {
  struct complex_t {
    cuFloatComplex v;
    __host__ __device__ complex_t( const cuFloatComplex& v ) : v(v) {}
    __host__ __device__ complex_t( const float_t& r, const float_t& i=0 ) : v(make_cuFloatComplex(r,i)) {}
    __host__ __device__ float real() const { return cuCrealf(v) ; } 
    __host__ __device__ float imag() const { return cuCimagf(v) ; } 
  };
  __host__ __device__ __inline__ complex_t operator+( complex_t a, complex_t b ){ return cuCaddf(a.v,b.v) ; } 
  __host__ __device__ __inline__ complex_t operator-( complex_t a, complex_t b ){ return cuCsubf(a.v,b.v) ; } 
  __host__ __device__ __inline__ complex_t operator/( complex_t a, complex_t b ){ return cuCdivf(a.v,b.v) ; } 
  __host__ __device__ __inline__ complex_t operator*( complex_t a, complex_t b ){ return cuCmulf(a.v,b.v) ; } 
  
  __host__ __device__ __inline__ complex_t operator+( float a, complex_t b ){ return complex_t( a + b.v.x, b.v.y ) ; }
  __host__ __device__ __inline__ complex_t operator-( float a, complex_t b ){ return complex_t( a - b.v.x, -b.v.y) ; } 
  __host__ __device__ __inline__ complex_t operator/( float a, complex_t b ){ float p = a / ( b.v.x*b.v.x + b.v.y*b.v.y) ; return complex_t( p * b.v.x, -p*b.v.y); } 
  __host__ __device__ __inline__ complex_t operator*( float a, complex_t b ){ return complex_t( a * b.v.x, a*b.v.y) ; } 
  
  __host__ __device__ __inline__ complex_t operator+( complex_t a, float b ){ return complex_t( a.v.x + b , a.v.y  ) ; } 
  __host__ __device__ __inline__ complex_t operator-( complex_t a, float b ){ return complex_t( a.v.x -b  , a.v.y ); } 
  __host__ __device__ __inline__ complex_t operator/( complex_t a, float b ){ return complex_t( a.v.x /b  , a.v.y/b ) ; } 
  __host__ __device__ __inline__ complex_t operator*( complex_t a, float b ){ return complex_t( b * a.v.x, b * a.v.y ) ; } 
  __host__ __inline__ std::ostream& operator<<( std::ostream& l, complex_t b ){ return l << "(" << b.real() << "," << b.imag() << ")"; }
};
using namespace ampgen_cuda; 

__global__ void bw_kernel( complex_t* __restrict__ r, const int N, const float* __restrict__ x0, const float3* __restrict__ x1)
{
  int i     = blockIdx.x * blockDim.x + threadIdx.x;
  float v899190686 = (sqrt((0.243717) + x1[i].x*x1[i].x + x1[i].y*x1[i].y + x1[i].z*x1[i].z) + sqrt((0.01948) + x1[i+N].x*x1[i+N].x + x1[i+N].y*x1[i+N].y + x1[i+N].z*x1[i+N].z))*(sqrt((0.243717) + x1[i].x*x1[i].x + x1[i].y*x1[i].y + x1[i].z*x1[i].z) + sqrt((0.01948) + x1[i+N].x*x1[i+N].x + x1[i+N].y*x1[i+N].y + x1[i+N].z*x1[i+N].z))-(x1[i].z + x1[i+N].z)*(x1[i].z + x1[i+N].z)-(x1[i].y + x1[i+N].y)*(x1[i].y + x1[i+N].y)-(x1[i].x + x1[i+N].x)*(x1[i].x + x1[i+N].x);
  float v2885272051 = (0.25)*v899190686-(0.131598) + (0.012571)/v899190686;
  r[i] = sqrt((0.900316)*x0[1]*x0[2]*x0[1]*sqrt(x0[2]*x0[2] + x0[1]*x0[1])/sqrt(x0[1]*sqrt(x0[2]*x0[2] + x0[1]*x0[1]) + x0[1]*x0[1]))*sqrt((9.)/((9.) + (3.)*v2885272051*x0[3]*x0[3] + v2885272051*x0[3]*x0[3]*v2885272051*x0[3]*x0[3]))/(x0[1]*x0[1]-v899190686-ampgen_cuda::complex_t(0.,1.)*x0[1]*x0[2]*((9.) + (3.)*x0[3]*((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1]))*x0[3] + x0[3]*((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1]))*x0[3]*x0[3]*((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1]))*x0[3])*sqrt(v2885272051/((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1])))*v2885272051*v2885272051*x0[1]*rsqrt(v899190686)/((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1]))/((0.25)*x0[1]*x0[1]-(0.131598) + (0.050282)/((4.)*x0[1]*x0[1]))/((9.) + (3.)*v2885272051*x0[3]*x0[3] + v2885272051*x0[3]*x0[3]*v2885272051*x0[3]*x0[3]));
}

#include "output.h"



std::vector<std::string> split( const std::string& s, char delim, bool ignoreWhitespace=true )
{
  std::vector<std::string> elems;
  std::string item;
  std::stringstream ss( s );
  while ( std::getline( ss, item, delim ) ) {
    if ( !ignoreWhitespace || ( item != " " && item != "" && item != "\n" && item != "\t" ) ) elems.push_back( item );
  }
  return elems;
}

int main(void)
{
  int N = 1 << 23; 
  complex_t * r;
  float     * pHost; 
  float3    * xE;
  std::ifstream stream("events.dat");
  std::string tmp;
  std::cout << sizeof(float3) << " " << sizeof(float) << std::endl; 
    
  cudaMallocManaged( &r    , sizeof(complex_t) * N  );
  cudaMallocManaged( &xE   , sizeof(float3)    * N * 3 );
  cudaMallocManaged( &pHost, sizeof(float)     * 4 );
  pHost[0] = 5 ;
  pHost[1] = 1.4324;
  pHost[2] = 0.109;
  pHost[3] = 1.5;

  std::vector<float> event_full( 12 * N );
  std::getline( stream , tmp );
  for (int i = 0; i < N+1; i++) {
    r[i]  = complex_t(0.0f,0.0f);
    std::getline( stream, tmp );
    auto tokens = split( tmp, ',');
      if( i== 0 ) std::cout << tmp << std::endl; 
    for( int p = 0 ; p < 3 ; ++p ){
      xE[ i + N * p  ].x  = stof( tokens[4*p]   );
      xE[ i + N * p  ].y  = stof( tokens[4*p+1] );
      xE[ i + N * p  ].z  = stof( tokens[4*p+2] );
    }
    for( int j = 0 ; j < 12; ++j ) 
      event_full[ j + 12*i ] = stof( tokens[j] );
  }
  for( int i = 0 ; i < 100 ; ++i )   bw_kernel<<< N/128, 128 >>>(r, N, pHost, xE );
  for( int i = 0 ; i < 100 ; ++i )   p2540321052 <<< N/128, 128 >>>(r, N, pHost, xE );

  cudaDeviceSynchronize();

  std::ofstream output("events_out.dat");
  for( int i = 0 ; i < N ; ++i )
  {
    for( int j = 0 ; j < 12 ; ++j ) output << event_full[j + 12*i] <<  " ";
    output << r[i] << std::endl;  
  }
  output.close();
  cudaFree(r);
  
  cudaFree(xE);
  cudaFree( pHost );
  return 0;
  
}
