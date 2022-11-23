#ifndef AMPGEN_SIMD_UTILS_H
#define AMPGEN_SIMD_UTILS_H

#include <array>
#include <complex>

#define INSTRUCTION_SET_SCALAR   0 
#define INSTRUCTION_SET_AVX2f    1 
#define INSTRUCTION_SET_AVX2d    2 
#define INSTRUCTION_SET_AVX512d  3 
#define INSTRUCTION_SET_ARM128d 10


#if INSTRUCTION_SET == INSTRUCTION_SET_SCALAR
  namespace scalar { using real_v = double; using complex_v = std::complex<double>; } 
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2f
//  #pragma message("Enable AVX2f")
  #include "AmpGen/simd/avx2f_types.h"
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2d
//  #pragma message("Enable AVX2d")
  #include "AmpGen/simd/avx2d_types.h"
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX512d
//  #pragma message("Enable AVX512d")
  #include "AmpGen/simd/avx512d_types.h"
#elif INSTRUCTION_SET == INSTRUCTION_SET_ARM128d
  #include "AmpGen/simd/arm128d_types.h"
#else 
  #pragma message("Unrecognised instruction set")
#endif

namespace AmpGen {

#if INSTRUCTION_SET == INSTRUCTION_SET_AVX512d
  namespace AVX        = AVX512d;
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2d
  namespace AVX        = AVX2d;
#elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2f
  namespace AVX        = AVX2f;
#elif INSTRUCTION_SET == INSTRUCTION_SET_SCALAR
  namespace AVX        = scalar; 
#elif INSTRUCTION_SET == INSTRUCTION_SET_ARM128d
  namespace AVX        = ARM128d; 
#endif

  using real_v      = AVX::real_v;
  using complex_v   = AVX::complex_v; 
  namespace utils {

    template <typename T> struct is_vector_type : std::false_type {}; 
    template <typename T> struct size           { static constexpr unsigned value = 1; } ;
    #if INSTRUCTION_SET != 0 
        template <>        struct is_vector_type <complex_v> : std::true_type {}; 
        template <>        struct is_vector_type <real_v  >  : std::true_type {}; 
        template <>        struct size           <complex_v>{ static constexpr unsigned value = real_v::size; }; 
        template <>        struct size           <real_v>   { static constexpr unsigned value = real_v::size; };
    
    #endif
    #if INSTRUCTION_SET == INSTRUCTION_SET_ARM128d
        template <>        struct size           <AVX::int_v>   { static constexpr unsigned value = 2; };
        template <>        struct is_vector_type <AVX::int_v>   : std::true_type {}; 
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
      if constexpr ( is_vector_type<simd_type>::value and std::is_same_v<simd_type, real_v> )
      {
        const auto arr = obj.to_ptr();
        auto rt = arr[0];
        for( unsigned i = 1 ; i != size<simd_type>::value; ++i ) rt += arr[i];
        return rt;
      }
      else return obj;
    }
    template <typename simd_type> bool all_of( const simd_type& obj)
    {
      if constexpr( size<simd_type>::value == 1 ) return obj;
      #if INSTRUCTION_SET == INSTRUCTION_SET_AVX2d
        return _mm256_movemask_pd( obj ) == 0xF;
      #elif INSTRUCTION_SET == INSTRUCTION_SET_AVX2f
        return _mm256_movemask_ps( obj ) == 0xFF;
      #endif
      return false; 
    }
    template <typename simd_type, typename value_type> bool all_of( const simd_type& obj, const value_type& v )
    {
      return all_of( obj == v );
    }
    template <typename T> auto make_complex( T&& re, T&& im ) { return std::complex<T>(re,im); }
    template <unsigned p=0, typename vtype> auto get( vtype v )
    { 
      if constexpr (  is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( std::is_same<vtype, complex_v>::value ) return std::complex( get<p>(v.real()), get<p>(v.imag()) );
      if constexpr ( !is_vector_type<vtype>::value ) return v; 
    }
    template < typename vtype> auto at( vtype v, const unsigned p=0 )
    {
      if constexpr ( is_vector_type<vtype>::value )
      {
        if constexpr ( std::is_same<vtype, real_v>::value ) return v.at(p);
        if constexpr ( std::is_same<vtype, complex_v>::value ) return std::complex( at( v.real(), p), at( v.imag(), p) );
      }
      else return v; 
    }
    template <typename>   struct is_std__complex : std::false_type{};
    template <typename T> struct is_std__complex<std::complex<T>> : std::true_type {} ; 

    template <typename T> inline auto norm( T&& value ){ 
      if constexpr( is_std__complex<std::remove_reference_t<T>> ::value ){ return std::norm(value); }
      else { return value.norm(); }
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
