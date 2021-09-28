#ifndef AMPGEN_KAHANSUM_H
#define AMPGEN_KAHANSUM_H

namespace AmpGen {
  /** @class KahanSum 
      @brief Implements Kahan summation for better precision with (repeated) floating point addition
  */
  template <typename T> 
  struct KahanSum {
    T sum = {0.};
    T cor = {0.};
    KahanSum& operator+=( const T& var )
    {
      T y = var - cor;
      T t = sum + y;
      cor = (t-sum)-y; 
      sum = t;
      return *this; 
    }
  };
}

#endif
