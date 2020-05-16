#ifndef AMPGEN_SIMD_UTILS_H
#define AMPGEN_SIMD_UTILS_H

#include <array>
#include <complex>

#if ENABLE_AVX512
  #include "AmpGen/simd/avx512d_types.h"
#elif ENABLE_AVX2d
  #include "AmpGen/simd/avx2d_types.h"
#elif ENABLE_AVX2f
  #include "AmpGen/simd/avx2f_types.h"
#endif



namespace AmpGen {
#if ENABLE_AVX512 
  namespace AVX        = AVX512d;  
#elif ENABLE_AVX2d
  namespace AVX        = AVX2d;
#elif ENABLE_AVX2f
  namespace AVX        = AVX2f;
#endif

#if ENABLE_AVX
  using float_v      = AVX::real_v;
  using complex_v    = AVX::complex_v;
#else 
  using float_v      = double;
  using complex_v    = std::complex<double>;
#endif

  namespace utils {

    template <class T> struct is_vector_type : std::false_type {}; 
    template <class T> struct size           { static constexpr unsigned value = 1; } ;
#if ENABLE_AVX
    template <>        struct is_vector_type <AVX::complex_v> : std::true_type {}; 
    template <>        struct is_vector_type <AVX::real_v  >  : std::true_type {}; 
    template <>        struct size           <AVX::complex_v>{ static constexpr unsigned value = AVX::complex_v::size; }; 
    template <>        struct size           <AVX::real_v>   { static constexpr unsigned value = AVX::complex_v::size; };
#endif
    template <class simd_type, class container_type, class functor_type> simd_type gather(
        const container_type& container, const functor_type& functor, unsigned offset=0, typename simd_type::scalar_type df =0.)
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
    template <class simd_type, class value_type> bool all_of( const simd_type& obj, const value_type& v )
    {
      if constexpr( ! is_vector_type<simd_type>::value ) return obj == v;
      else return _mm256_movemask_pd( obj == v ) == 0xF;
    }
    template <class simd_type> bool all_of( const simd_type& obj)
    {
      if constexpr( ! is_vector_type<simd_type>::value ) return obj;
      else return _mm256_movemask_pd( obj ) == 0xF;
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
      #if ENABLE_AVX 
      if constexpr(   is_vector_type<ctype>::value ) return AVX::norm(v);
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
