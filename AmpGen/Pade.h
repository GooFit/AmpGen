#ifndef AMPGEN_PADE_H
#define AMPGEN_PADE_H

#include <array>
#include <functional>
#include <vector>
#include <iostream>

namespace AmpGen { 
  enum Strategy { linear, quadratic, cubic, quartic };
  namespace detail { 
    std::vector<double> solve_pade(const std::function<double(const double&)>& fcn,
        const double& min,
        const double& max,
        const unsigned& N,
        const Strategy& strat = Strategy::linear);
  }

  template <unsigned N, class T = double> class Pade 
  {
    public:
      Pade( const std::vector<double>& r, 
            const double& min, 
            const double& max) : min(min), max(max) 
      {
        for(unsigned i = 0; i <= N; ++i ) co_f[i] = r[i];
        for(unsigned i = 0; i <  N; ++i ) co_g[i] = r[i+(N+1)];
        range = 1./(max-min);
      }
      
      Pade(const std::function<double(const double&)>& fcn, 
          const double& min, 
          const double& max,
          const Strategy& strat = Strategy::linear) : 
        m_function(fcn), min(min),max(max)
    {
      auto r = detail::solve_pade(fcn, min, max, N, strat );
      for(unsigned i = 0; i <= N; ++i ) co_f[i] = r[i];
      for(unsigned i = 0; i <  N; ++i ) co_g[i] = r[i+(N+1)];
      range = 1./(max-min);
    }
      template <typename T2> 
      T2 operator()(const T2& s) const 
      {
        T2 x = (s-min)*range;
        T2 f = 0.;
        T2 g = 1.;
        T2 acc = 1.;
        for(unsigned i = 0; i < N; ++i){
          f += co_f[i] * acc;
          acc *= x;
          g += co_g[i] * acc;
        }
        return (f + co_f[N]*acc)/g;
      }
      void print() const {
        for( int i = 0 ; i != N+1; ++i ) std::cout << co_f[i] << std::endl; 
        for( int i = 0 ; i != N ; ++i ) std::cout << co_g[i] << std::endl; 
      }
    private:
      std::function<double(const double&)> m_function;
      std::array<T, N+1> co_f;
      std::array<T, N>   co_g; 
      T min;
      T max;
      T range; 
  };
}

#endif
