#include "AmpGen/Pade.h"
#include "TMatrixD.h"
#include "TVectorD.h"


std::vector<double> AmpGen::detail::solve_pade(const std::function<double(const double&)>& fcn,
    const double& min,
    const double& max,
    const unsigned& N, 
    const AmpGen::Strategy& strat)
{

  TMatrixD solver(2*N+1,2*N+1);
  std::vector<double> samples(2*N+1);
  if( strat < 4 ){
    for(unsigned eq = 0 ; eq < 2*N+1; ++eq) 
      samples[eq] = pow( eq/double(2*N), strat + 1);
  }
  TVectorD rest(2*N+1); 
  for( unsigned eq = 0 ; eq < 2*N+1; ++eq ){
    rest(eq) = fcn( samples[eq]*(max-min) + min);
    for(unsigned i = 0; i <= N; ++i) solver(eq,i)   = pow(samples[eq],i);
    for(unsigned i = 1; i <= N; ++i) solver(eq,i+N) = -rest(eq)* pow(samples[eq],i);
  }
  solver.Invert();
  auto r = solver * rest; 
  std::vector<double> rv(2*N+1);
  for(unsigned i = 0; i<2*N+1; ++i) rv[i] = r[i];
  return rv;
}
