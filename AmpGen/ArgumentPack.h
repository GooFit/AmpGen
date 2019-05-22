#ifndef AMPGEN_ARGUMENTPACK_H
#define AMPGEN_ARGUMENTPACK_H

#include "AmpGen/MetaUtils.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace AmpGen
{
  struct IArgument {
    virtual ~IArgument() = default;
    
  };
  /** @class Argument 
    @brief Structure to pass "named" parameters to functions. 
    Structure to flexibly pass blocks of "Named" parameters to functions, to 
    approximate the behaviour of Python's named arguments. 
    Typical usage is for constructors with variable arguments, such as 
    to read data from the disk. The interface for the user is typically 
    \code{cpp}
    EventList events("files.root", 
    type,
    GetGenPDF(true), 
    WeightBranch("eventWeight"),
    Branches({"K_PX","K_PY",...}));
    \endcode  
    Internally these arguments are used to construct an ArgumentPack, and then read out 
    using getArg<TYPE>, so for this example: 
    \code{cpp}
    auto pdfSize      = args.getArg<CacheSize>().val;
    auto filter       = args.getArg<Filter>().val;
    auto getGenPdf    = args.getArg<GetGenPdf>(true).val;
    auto weightBranch = args.getArg<WeightBranch>().val;
    auto branches     = args.getArg<Branches>().val;
    auto applySym     = args.getArg<ApplySym>().val;
    auto entryList    = args.getArg<EntryList>().val; 
    \endcode
    @tparam TYPE Type of the argument, such as a string, a number, a bool etc.  
  */
  template <typename TYPE>
  struct Argument : public IArgument 
  {
    template <typename T>
    Argument( T x ) : val(x) {}
    Argument() {
      if constexpr(std::is_pointer_v<TYPE>) val = nullptr;
      else val = TYPE();
    }
    operator TYPE() const { return val; }
    TYPE val;
  };
  
  struct File : public IArgument 
  {
    std::string name;
    std::ios_base::openmode mode;
    File( const std::string& name = "", const std::ios_base::openmode& mode = std::ios_base::in )
      : name( name ), mode( mode ){};
  };

  class ArgumentPack
  {
    public:
      template <typename... ARGS>
        ArgumentPack( const ARGS&... args )
        {
          std::tuple<ARGS...> argTuple( args... );
          for_each( argTuple, [this]( auto& f ) { m_parameters.emplace_back( makeShared( f ) ); } );
        }
      template <typename ARG, typename DEFAULT_TYPE=ARG> 
        ARG getArg( const DEFAULT_TYPE& default_argument = DEFAULT_TYPE() ) const
      {
        for ( auto& param : m_parameters ) {
          auto ptr = dynamic_cast<ARG*>( param.get() );
          if ( ptr != nullptr ) return *ptr;
        }
        return ARG(default_argument);
      }
    private:
      std::vector<std::shared_ptr<IArgument>> m_parameters;
  };
#define DECLARE_ARGUMENT(X, Y)                            \
  struct X : public AmpGen::Argument<Y> {                 \
    template<class Z>                                     \
    X(Z val) : AmpGen::Argument<Y>(val){}                 \
    X() : AmpGen::Argument<Y>(){}                         \
  }
} // namespace AmpGen

#endif
