#ifndef AMPGEN_LITESPAN_H
#define AMPGEN_LITESPAN_H

#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

namespace AmpGen {
  // replace with std::span when c++20 becomes widespread 
  template <typename type> class LiteSpan 
  {
    public:
      LiteSpan( const type* data, unsigned size =1) : m_data(data), m_size(size) {}    
      const type& operator[](const unsigned index) const { return m_data[index]; }
      operator type() const { return m_data[0] ; }
      unsigned size() const { return m_size; }
    private:
      const type* m_data = {nullptr};
      unsigned    m_size = {0}; 
  };

  template <typename return_type, typename arg_type> struct KeyedFunctors 
  {
    std::vector<std::function<return_type(const arg_type&)> > functors; 
    std::vector<std::string> keys;
    std::vector<std::string> titles;
    template <class functor_type> 
    void add(const functor_type& functor, const std::string& key, const std::string& title="") 
    {
      functors.push_back(functor);
      keys.push_back(key);
      titles.push_back(title);
    }
    std::vector<return_type> operator()( const arg_type& arg ) const { 
      
      std::vector<return_type> rt; 
      for( auto& f : functors ) rt.push_back( f(arg) );
      return rt; }
  };

}

#endif
