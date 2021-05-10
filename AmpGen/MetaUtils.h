#ifndef AMPGEN_METAUTILS_H
#define AMPGEN_METAUTILS_H

#include <memory>
#include <tuple>
#include <vector>
#include <string>
#include <cxxabi.h>
#include <algorithm>

namespace AmpGen
{
  /** Utility classes for compile-time metaprogramming, such as identifying the types of 
    arguments for generating source code, compile-time unrolling of tuples and loops, 
    and identifying if a class can be constructed in different ways. 
    */

  template <class TYPE> std::string type_string()
  {
    int status = 0;
    std::string name = abi::__cxa_demangle( typeid( TYPE ).name(), nullptr, nullptr, &status );
    if( std::is_const<typename std::remove_reference<TYPE>::type>::value ) 
      name = "const " + name;
    if( std::is_reference<TYPE>() ) name = name + "&";
    return name; 
  }

  template <class TYPE> std::string type_string( const TYPE& t ) { return type_string<TYPE>(); }

  namespace detail {
    template<typename T, typename... args> struct zeroType { typedef T type; };
  }
  template<int N, typename... args> using nthType = typename std::tuple_element<N, std::tuple<args...>>::type;

  template<typename... args> using zeroType = typename detail::zeroType<args...>::type; 

  template <std::size_t I = 0, typename FuncT, typename... Tp>
    typename std::enable_if_t<I == sizeof...( Tp ), void>
    for_each( std::tuple<Tp...>&, FuncT ){}

  template <std::size_t I = 0, typename FuncT, typename... Tp> 
    inline typename std::enable_if_t< I<sizeof...( Tp ), void> for_each( std::tuple<Tp...>& t, FuncT f )
    {
      f( std::get<I>( t ) );
      for_each<I + 1, FuncT, Tp...>( t, f );
    }
  template <typename iterator, typename... transform_types> 
    void for_each_sequence( iterator begin, iterator end, transform_types... transforms)
    {
      for_each( std::tuple<transform_types...>(transforms...), [&](auto& transform){ std::for_each( begin, end, transform ); } );
    }

  template <std::size_t I = 0, typename FuncT, typename... Tp>
    typename std::enable_if_t<I == sizeof...( Tp ), void>
    for_each( const std::tuple<Tp...>&, FuncT ){}

  template <std::size_t I = 0, typename FuncT, typename... Tp> 
    inline typename std::enable_if_t< I<sizeof...( Tp ), void> for_each( const std::tuple<Tp...>& t, FuncT f )
    {
      f( std::get<I>( t ) );
      for_each<I + 1, FuncT, Tp...>( t, f );
    }

  template <typename R, 
           typename RT, 
           typename ARGS, 
           bool result = std::is_same<decltype( ( ( *(R*)nullptr ) )( ARGS() ) ), RT>::value>
             constexpr bool typeHelper( int )
             {
               return result;
             }

  template <typename R, 
            typename RT, 
            typename ARGS> 
             constexpr bool typeHelper( ... )
             {
               return false;
             }

  template <typename R, 
           typename RT, 
           typename ARGS> 
             constexpr bool hasReturnType()
             {
               return typeHelper<R, RT, ARGS>( 0 );
             }
  template <typename T, typename... R> constexpr bool hasConstructor()
  {
    return std::is_constructible<T,R...>::value && (false == std::is_same<T, R...>::value);
  }


  template <typename arg=void, typename... args> std::vector<std::string> typelist()
  {
    std::vector< std::string > rt;
    if( type_string<arg>() != "void" ) {
      rt.emplace_back( type_string<arg>() );
      auto rtp = typelist<args...>(); 
      std::copy( rtp.begin(), rtp.end(), std::back_inserter(rt) );
    }
    return rt;
  }

  template <typename> struct isTuple: std::false_type {};
  template <typename ...T> struct isTuple<std::tuple<T...>>: std::true_type {};
  template <typename> struct isVector : std::false_type {};
  template <typename    T> struct isVector<std::vector<T>> : std::true_type {};
  template <typename> struct is_array : std::false_type {};
  template <typename    T, std::size_t N> struct is_array<std::array<T,N>> : std::true_type {};


#define def_has_function(function_name) \
  template <typename T>      \
  struct has_##function_name { \
    template<typename U> static auto test(int) -> decltype(std::declval<U>().function_name() == 1, std::true_type()); \
    template<typename>   static std::false_type test(...); \
    static constexpr bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value; \
  };

  template<typename ...T>
    struct is_functor : std::false_type {};

  template<typename T, typename rt, typename... args>
    struct is_functor<T, rt(args...)> {
        template<typename U>
        static constexpr auto test(U*) -> decltype(std::declval<U>().operator()( std::declval<args>()...), std::true_type()); 
        template<typename> static constexpr std::false_type test(...);
        static constexpr bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value ;
    };

  template <int A, int B> struct get_power
  {
    static const int value = A * get_power<A, B - 1>::value;
  };
  template <int A> struct get_power<A, 0>
  {
    static const int value = 1;
  };

} // namespace AmpGen

#endif
