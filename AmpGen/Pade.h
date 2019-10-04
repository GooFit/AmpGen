#ifndef AMPGEN_PADE_H
#define AMPGEN_PADE_H

#include <array>
#include <functional>
#include <vector>

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
      T operator()(const T& s)
      {
        T x = (s-min)*range;
        T f = 0;
        T g = 1;
        T acc = 1;
        for(unsigned i = 0; i < N; ++i){
          f += co_f[i] * acc;
          acc *= x;
          g += co_g[i] * acc;
        }
        return (f + co_f[N]*acc)/g;
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
