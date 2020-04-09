#include <array>
#include "AmpGen/simd/avx2_types.h"

namespace AmpGen {
  namespace utils {

    template <class T> struct is_vector_type { static constexpr bool value = false; }; 
    template <>        struct is_vector_type <AVX2::complex_t>{ static constexpr bool value = true ; };
    template <>        struct is_vector_type <AVX2::float_t  >{ static constexpr bool value = true ; };

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
      return simd_type::size * unsigned ( 1 + (unaligned_size -1 ) / simd_type::size );
    }
    template <class simd_type> auto sum_elements( const simd_type& obj )
    {
      if constexpr ( is_vector_type<simd_type>::value )
      {
        auto arr = obj.to_array();
        auto rt = arr[0];
        for( unsigned i = 1 ; i != simd_type::size; ++i ) rt = rt + arr[i];
        return rt;
      }
      else return obj;
    }
    template <unsigned p=0, class vtype> auto get( vtype v ){ 
      if constexpr ( is_vector_type<vtype>::value ) return v.at(p); 
      if constexpr ( ! is_vector_type<vtype>::value ) return v; 
    } 

  }
}
