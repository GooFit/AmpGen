#ifndef AMPGEN_SIMD_UTILS_H
#define AMPGEN_SIMD_UTILS_H

#include <array>
#if ENABLE_AVX2 
#if DOUBLE_PRECISION
  #include "AmpGen/simd/avx2d_types.h"
#else
  #include "AmpGen/simd/avx2_types.h"
#endif
#endif

namespace AmpGen {
#if ENABLE_AVX2
#if DOUBLE_PRECISION
  using float_v    = AVX2d::float_t;
  using complex_v  = AVX2d::complex_t;
#else
  using float_v    = AVX2::float_t;
  using complex_v  = AVX2::complex_t;
#endif
#else 
  using float_v    = double;
  using complex_v  = std::complex<double>;
#endif

  namespace utils {

    template <class T> struct is_vector_type : std::false_type {}; 
    template <class T> struct size           { static constexpr unsigned value = 1; } ;
#if ENABLE_AVX2 
#if DOUBLE_PRECISION
    template <>        struct is_vector_type <AVX2d::complex_t> : std::true_type {}; 
    template <>        struct is_vector_type <AVX2d::float_t  > : std::true_type {}; 
    template <>        struct size           <AVX2d::complex_t>{ static constexpr unsigned value = 4; }; 
    template <>        struct size           <AVX2d::float_t>  { static constexpr unsigned value = 4; };
#else
    template <>        struct is_vector_type <AVX2::complex_t> : std::true_type {}; 
    template <>        struct is_vector_type <AVX2::float_t  > : std::true_type {}; 
    template <>        struct size           <AVX2::complex_t>{ static constexpr unsigned value = 8; }; 
    template <>        struct size           <AVX2::float_t>  { static constexpr unsigned value = 8; };
#endif
#else
    template <>        struct is_vector_type <float_v>   : std::false_type {}; 
    template <>        struct is_vector_type <complex_v> : std::false_type {}; 
#endif
    template <class simd_type, class container_type, class functor_type> simd_type gather(
        const container_type& container, const functor_type& functor, unsigned offset=0, float df =0.)
    {
      std::array<typename simd_type::scalar_type, simd_type::size> rv; 
      if( df == 0. )
      for( unsigned k = 0 ; k != simd_type::size; ++k ) rv[k] = offset + k < container.size() ? functor(container[offset+k]) : functor(container[container.size()-1]); 
      else 
      for( unsigned k = 0 ; k != simd_type::size; ++k ) rv[k] = offset + k < container.size() ? functor(container[offset+k]) : df; 
      return simd_type( rv.data() );
    } 

    template <class simd_type> size_t aligned_size( const size_t& unaligned_size ) {
      return size<simd_type>::value * unsigned ( 1 + (unaligned_size -1 ) / size<simd_type>::value );
    }
    template <class simd_type> auto sum_elements( const simd_type& obj )
    {
      if constexpr ( is_vector_type<simd_type>::value )
      {
        auto arr = obj.to_array();
        auto rt = arr[0];
        for( unsigned i = 1 ; i != size<simd_type>::value; ++i ) rt = rt + arr[i];
        return rt;
      }
      else return obj;
    }
    template <unsigned p=0, class vtype> auto get( vtype v )
    { 
      if constexpr ( is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( ! is_vector_type<vtype>::value ) return v; 
    }
    template < class vtype> auto at( vtype v, const unsigned p=0 )
    {
      if constexpr ( is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( ! is_vector_type<vtype>::value ) return v; 
    }
    template <class ctype> auto norm( const ctype& v )
    {
      #if ENABLE_AVX2 && DOUBLE_PRECISION 
      if constexpr(   is_vector_type<ctype>::value ) return AVX2d::norm(v);
      #endif 
      #if ENABLE_AVX2 && ! DOUBLE_PRECISION 
      if constexpr(   is_vector_type<ctype>::value ) return AVX2::norm(v);
      #endif 
      if constexpr( ! is_vector_type<ctype>::value ) return std::norm(v);
    }
    template <class type, class store_type> void store( store_type* container, const type& v)
    {
      if constexpr( is_vector_type<type>::value ) 
      {
        auto arr = v.to_array();
        for( unsigned k = 0 ; k != utils::size<type>::value; ++k ) container[k] =  arr[k];
      }
      else {
        *container = v;
      }
    }
  }
}

#endif
