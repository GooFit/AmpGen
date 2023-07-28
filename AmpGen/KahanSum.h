#ifndef AMPGEN_KAHANSUM_H
#define AMPGEN_KAHANSUM_H

#ifdef __FAST_MATH__
#error Compensated summation is unsafe with -ffast-math (/fp:fast)
#endif
#include "AmpGen/simd/utils.h"

namespace AmpGen {
  /** @class KahanSum 
    @brief Implements Kahan summation for better precision with (repeated) floating point addition
    */
  template <typename T> 
    struct KahanSum {
      T sum = 0.;
      T cor = 0.;
      KahanSum& operator+=( const T& var )
      {
        
        T t = sum + var;
        if constexpr( utils::is_vector_type<T>::value )
        { 
          cor += select( abs(sum) >= abs(var), (sum-t)+var, (var-t) + sum ); 
        }
        else {
         
          cor += std::abs(sum) >= std::abs(var) ?  (sum-t)+var :  (var-t) + sum; 
        }
        sum = t;
        return *this; 
      }
      T operator()() const { return sum + cor; } 
    };
  template <typename T> 
    KahanSum<T> operator+( const KahanSum<T>& l, const T& var )
    {
      KahanSum<T> rt; 
      T y = var - l.cor;
      T t = l.sum + y;
      rt.cor = (t-l.sum)-y; 
      rt.sum = t;
      return rt;
    }


  template <typename T> T KahanBabushkaKleinSum(const std::vector<T>& container)
  {
    T sum = 0.0;
    T cs  = 0.0;
    T ccs = 0.0;
    T c   = 0.0;
    T cc  = 0.0;

    for( auto& input : container )
    {
      T t = sum + input;
      if ( std::fabs(sum) >= std::fabs(input) ) 
      {
        c = (sum - t) + input;
      }
      else c = (input - t) + sum;
      sum = t;
      t = cs + c;
      if ( std::fabs(cs) >= std::fabs(c) )
        cc = (cs - t) + c;
      else
        cc = (c - t) + cs;
      cs = t;
      ccs = ccs + cc;
    }
    return sum + cs + ccs;
  }
}

#endif
