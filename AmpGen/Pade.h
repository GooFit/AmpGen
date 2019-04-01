#ifndef AMPGEN_PADE_H
#define AMPGEN_PADE_H

#include "TMatrixD.h"
#include "TVectorD.h"
#include <array>
#include <functional>

namespace AmpGen { 
  enum Strategy { linear, quadratic, cubic, quartic };
  template <int N, class T = double>
    class Pade {
      private:
        std::function<double(const double&)> m_function;
        std::array<T, N+1> co_f;
        std::array<T, N>   co_g; 
        T min;
        T max;
        T range; 
      public:

        Pade(const std::function<double(const double&)>& fcn, 
            const double& min, 
            const double& max,
            const Strategy& strat = Strategy::linear) : 
          m_function(fcn), min(min),max(max){
            TMatrixD solver(2*N+1,2*N+1);
            std::vector<double> samples(2*N+1);
            if( strat < 4 ){
              for(size_t eq = 0 ; eq < 2*N+1; ++eq) 
                samples[eq] = pow( eq/double(2*N), strat + 1);
            }
            TVectorD rest(2*N+1); 
            for( int eq = 0 ; eq < 2*N+1; ++eq ){
              rest(eq) = fcn( samples[eq]*(max-min) + min);
              for(int i = 0; i <= N; ++i) solver(eq,i)   = pow(samples[eq],i);
              for(int i = 1; i <= N; ++i) solver(eq,i+N) = -rest(eq)* pow(samples[eq],i);
            }
            solver.Invert();
            auto r = solver * rest; 
            for(size_t i = 0; i <= N; ++i ) co_f[i] = r[i];
            for(size_t i = 0; i <  N; ++i ) co_g[i] = r[i+(N+1)];
            range = 1./(max-min);
          }
        T operator()(const T& s){
          T x = (s-min)*range;
          T f = 0;
          T g = 1;
          T acc = 1;
          for(size_t i = 0; i < N; ++i){
            f += co_f[i] * acc;
            acc *= x;
            g += co_g[i] * acc;
          }
          return (f + co_f[N]*acc)/g;
        }
    };
}

#endif
