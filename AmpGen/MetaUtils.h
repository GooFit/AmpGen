#ifndef AMPGEN_METAUTILS_H
#define AMPGEN_METAUTILS_H

#include <memory>
#include <tuple>
#include <vector>
#include <string>
#include <cxxabi.h>

namespace AmpGen
{
  template <class TYPE> std::string typeof()
  {
    int status = 0;
    std::string name = abi::__cxa_demangle( typeid( TYPE ).name(), nullptr, nullptr, &status );
    if( std::is_const<typename std::remove_reference<TYPE>::type>::value ) 
      name = "const " + name;
    if( std::is_reference<TYPE>() ) name = name + "&";
    return name; 
  }

  template <class TYPE> std::string typeof( TYPE t ) { return typeof<TYPE>(); }

  template <std::size_t I = 0, typename FuncT, typename... Tp>
  typename std::enable_if<I == sizeof...( Tp ), void>::type
  for_each( std::tuple<Tp...>&, FuncT ){}

  template <std::size_t I = 0, typename FuncT, typename... Tp> 
  inline typename std::enable_if < I<sizeof...( Tp ), void>::type for_each( std::tuple<Tp...>& t, FuncT f )
  {
    f( std::get<I>( t ) );
    for_each<I + 1, FuncT, Tp...>( t, f );
  }

  template <class TYPE> std::shared_ptr<TYPE> makeShared( const TYPE& obj )
  {
    return std::make_shared<TYPE>( obj );
  }

  template <typename R, typename RT, typename ARGS, bool result = std::is_same<decltype( ( ( *(R*)nullptr ) )( ARGS() ) ), RT>::value>
  constexpr bool typeHelper( int )
  {
    return result;
  }

  template <typename R, typename RT, typename ARGS> constexpr bool typeHelper( ... )
  {
    return false;
  }

  template <typename R, typename RT, typename ARGS> constexpr bool hasReturnType()
  {
    return typeHelper<R, RT, ARGS>( 0 );
  }

  template <typename arg=void, typename... args> std::vector<std::string> typelist()
  {
    std::vector< std::string > rt;
    if( typeof<arg>() != "void" ) {
      rt.emplace_back( typeof<arg>() );
      auto rtp = typelist<args...>(); 
      for( auto& r : rtp ) rt.emplace_back( r );
    }
    return rt;
  }
} // namespace AmpGen

#endif
