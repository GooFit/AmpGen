#ifndef AMPGEN_LITESPAN_H
#define AMPGEN_LITESPAN_H

#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

namespace AmpGen {
  template <typename return_type, typename ...arg_types> struct KeyedFunctors;

  template <typename return_type, typename ...arg_types> struct KeyedFunctors <return_type(arg_types...)>
  {
    std::vector<std::function<return_type(arg_types...)> > functors; 
    std::vector<std::string> keys;
    std::vector<std::string> titles;
    template <typename functor_type> void add(const functor_type& functor, const std::string& key, const std::string& title="") 
    {
      functors.push_back(functor);
      keys.push_back(key);
      titles.push_back(title);
    }
    std::vector<return_type> operator()( arg_types... arg ) const 
    {  
      std::vector<return_type> rt; 
      for( auto& f : functors ) rt.push_back( f(arg...) );
      return rt;
    }
  };

}

#endif
