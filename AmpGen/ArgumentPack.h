#ifndef AMPGEN_ARGUMENTPACK_H
#define AMPGEN_ARGUMENTPACK_H

#include "AmpGen/MetaUtils.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>


#define DECLARE_ARGUMENT( X, Y )                                                                                       \
  struct X : public AmpGen::Argument<Y> {                                                                              \
    X( const Y& val = Y() ) : AmpGen::Argument<Y>( val ){};                                                            \
  }

#define DECLARE_ARGUMENT_PTR( X, Y )                                                                                   \
  struct X : public AmpGen::ArgumentPtr<Y> {                                                                           \
    X( Y* val = nullptr ) : AmpGen::ArgumentPtr<Y>( val ){};                                                           \
  }

#define DECLARE_ARGUMENT_DEFAULT( X, Y, Z )                                                                            \
  struct X : public AmpGen::Argument<Y> {                                                                              \
    X( const Y& val = Z ) : AmpGen::Argument<Y>( val ){};                                                              \
  }

namespace AmpGen
{

  struct IArgument {
    virtual ~IArgument() = default;
  };

  template <class TYPE>
  struct Argument : public IArgument {
    Argument( const TYPE& x = TYPE() ) : val( x ) {}
    TYPE val;
    operator TYPE() const { return val; }
  };
  template <class TYPE>
  struct ArgumentPtr : public IArgument {
    ArgumentPtr( TYPE* x = nullptr ) : val( x ) {}
    TYPE* val;
  };
  struct File : public IArgument {
    std::string name;
    std::ios_base::openmode mode;
    File( const std::string& name = "", const std::ios_base::openmode& mode = std::ios_base::in )
        : name( name ), mode( mode ){};
  };

  class ArgumentPack
  {

  public:
    template <class ARG>
    ARG getArg( const ARG& default_argument = ARG() ) const
    {
      for ( auto& param : m_parameters ) {
        auto ptr = dynamic_cast<ARG*>( param.get() );
        if ( ptr != nullptr ) return *ptr;
      }
      return default_argument;
    }

    template <class... ARGS>
    ArgumentPack( const ARGS&... args )
    {
      std::tuple<ARGS...> argTuple( args... );
      for_each( argTuple, [this]( auto& f ) { m_parameters.emplace_back( makeShared( f ) ); } );
    }
    std::vector<std::shared_ptr<IArgument>> m_parameters;
  };
} // namespace AmpGen

#endif
