#ifndef AMPGEN_SIMD_UTILS_H
#define AMPGEN_SIMD_UTILS_H

#include <array>
#include <complex>

namespace AmpGen {
  namespace AVX2d  { class real_v; class complex_v; }
  namespace AVX2f  { class real_v; class complex_v; }
  namespace AVX512 { class real_v; class complex_v; }
}

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

    template <typename T> struct is_vector_type : std::false_type {}; 
    template <typename T> struct size           { static constexpr unsigned value = 1; } ;
#if ENABLE_AVX
    template <>        struct is_vector_type <AVX::complex_v> : std::true_type {}; 
    template <>        struct is_vector_type <AVX::real_v  >  : std::true_type {}; 
    template <>        struct size           <AVX::complex_v>{ static constexpr unsigned value = AVX::complex_v::size; }; 
    template <>        struct size           <AVX::real_v>   { static constexpr unsigned value = AVX::complex_v::size; };
#endif
    template <typename simd_type, typename container_type, typename functor_type> simd_type gather(
        const container_type& container, const functor_type& functor, unsigned offset=0, typename simd_type::scalar_type df =0.)
    {
      std::array<typename simd_type::scalar_type, simd_type::size> rv; 
      if( df == 0. )
      for( unsigned k = 0 ; k != simd_type::size; ++k ) rv[k] = offset + k < container.size() ? functor(container[offset+k]) : functor(container[container.size()-1]); 
      else 
      for( unsigned k = 0 ; k != simd_type::size; ++k ) rv[k] = offset + k < container.size() ? functor(container[offset+k]) : df; 
      return simd_type( rv.data() );
    } 

    template <typename simd_type> size_t aligned_size( const size_t& unaligned_size ) {
      return size<simd_type>::value * unsigned ( 1 + (unaligned_size -1 ) / size<simd_type>::value );
    }
    template <typename simd_type> auto sum_elements( const simd_type& obj )
    {
      if constexpr ( is_vector_type<simd_type>::value )
      {
        if constexpr( std::is_same_v<simd_type, AVX2d::real_v> or std::is_same_v<simd_type, AVX2f::real_v> )
        {
          const auto arr = obj.to_ptr();
          auto rt = arr[0];
          for( unsigned i = 1 ; i != size<simd_type>::value; ++i ) rt = rt + arr[i];
          return rt;
        }
        if constexpr( std::is_same_v<simd_type, AVX2d::complex_v> or std::is_same_v<simd_type, AVX2f::complex_v> )
        {
          const auto arr = obj.to_array();
          auto rt = arr[0];
          for( unsigned i = 1 ; i != size<simd_type>::value; ++i ) rt = rt + arr[i];
          return rt;
        }
      }
      else return obj;
    }
    template <typename simd_type> bool all_of( const simd_type& obj)
    {
      if constexpr( size<simd_type>::value == 1 ) return obj;
      if constexpr( std::is_same_v<simd_type, AVX2d::real_v> ) return _mm256_movemask_pd( obj ) == 0xF;
      if constexpr( std::is_same_v<simd_type, AVX2f::real_v> ) return _mm256_movemask_ps( obj ) == 0xFF;
      return false; 
    }
    template <typename simd_type, typename value_type> bool all_of( const simd_type& obj, const value_type& v )
    {
      return all_of( obj == v );
    }
    template <unsigned p=0, typename vtype> auto get( vtype v )
    { 
      if constexpr ( is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( ! is_vector_type<vtype>::value ) return v; 
    }
    template < typename vtype> auto at( vtype v, const unsigned p=0 )
    {
      if constexpr ( is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( ! is_vector_type<vtype>::value ) return v; 
    }
    template <typename ctype> auto norm( const ctype& v )
    {
      #if ENABLE_AVX 
      if constexpr(   is_vector_type<ctype>::value ) return AVX::norm(v);
      #endif 
      if constexpr( ! is_vector_type<ctype>::value ) return std::norm(v);
    }
    template <typename type, typename store_type> void store( store_type* container, const type& v)
    {
      if constexpr( is_vector_type<type>::value ) 
      {
        auto arr = v.to_ptr();
        for( unsigned k = 0 ; k != utils::size<type>::value; ++k ) container[k] =  arr[k];
      }
      else {
        *container = v;
      }
    }
  }
}

#endif
