#include <cmath>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "AmpGen/Pade.h"

std::vector<double> AmpGen::detail::solve_pade(const std::function<double(const double&)>& fcn,
    const double& min,
    const double& max,
    const unsigned& N, 
    const AmpGen::Strategy& strat)
{
  gsl_matrix* solver  = gsl_matrix_alloc( 2*N+1, 2*N+1);
  gsl_matrix* inverse = gsl_matrix_alloc( 2*N+1, 2*N+1);
  gsl_matrix_set_all(solver, 0. );

  std::vector<double> samples(2*N+1, 0);
  std::vector<double> rest( 2*N+1, 0);
  for(unsigned eq = 0 ; eq < 2*N+1; ++eq) 
    samples[eq] = std::pow( eq/double(2*N), strat + 1);
  for( unsigned eq = 0 ; eq < 2*N+1; ++eq ){
    rest[eq] = fcn( samples[eq]*(max-min) + min);
    for(unsigned i = 0; i <= N; ++i) gsl_matrix_set(solver, eq, i, std::pow(samples[eq],i) );
    for(unsigned i = 1; i <= N; ++i) gsl_matrix_set(solver, eq, i+N, -rest[eq]* std::pow(samples[eq],i) );
  }
  gsl_permutation *p = gsl_permutation_alloc(2*N+1);
  int s;
  gsl_linalg_LU_decomp(solver, p, &s);
  gsl_linalg_LU_invert(solver, p, inverse); 
  
  std::vector<double> rv(2*N+1,0);
  for(unsigned i = 0; i<2*N+1; ++i){ 
    for( unsigned j = 0 ; j != 2*N+1; ++j ) rv[i] += gsl_matrix_get(inverse, i,j) * rest[j];
  }
  gsl_permutation_free(p);
  gsl_matrix_free(solver); 
  gsl_matrix_free(inverse); 
/*  
  std::vector<double> co_f( N + 1);
  std::vector<double> co_g( N ); 
  for( int i = 0 ; i != N; ++i )
  {
    co_f[i] = rv[i];
    co_g[i] = rv[i+(N+1)];
  }
  co_f[N] = rv[N];

  auto test_fcn = [co_f, co_g, max, min, N ]( double s )
  {
    auto range = 1/(max-min);
    auto x = ( s -min ) * range; 
    auto f = 0.;
    auto g = 1.;
    auto acc = 1.;
    for(unsigned i = 0; i < N; ++i){
      f += co_f[i] * acc;
      acc *= x;
      g += co_g[i] * acc;
    }
    return (f + co_f[N]*acc)/g;
  };
  for( int j = 0 ; j != 2*N+1; ++j )
  {
    double s = samples[j] * ( max - min ) + min;
    std::cout << j << " " << s << " " << test_fcn(s) - fcn(s) << std::endl; 
  }
  */
  return rv;
}
