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
  /// functor-like object that documents what is stored in each slot;
  /// This mutated into a cache-like object as was I was writing it, so 
  /// should rename it to something else ...  
  template <typename return_type, 
            typename container_type,
            typename cache_type = return_type> class KeyedView
  {
      typedef typename container_type::value_type value_type;
    public:
      KeyedView( const container_type& container, const unsigned width ) : 
        m_container(&container), 
        m_cache( width * container.size(),0 ), 
        m_width(width),
        m_size(container.size()),
        m_keys( width, "") {}
      unsigned index(const value_type& it)                                      const { return it.index() ; }
      const std::string& key(const unsigned int& column )                       const { return m_keys[column] ; }
      const return_type* operator()( const value_type& it )                     const { 
        if( m_width *index(it) >= m_cache.size()) ERROR("Out-of-bounds access : " << index(it) );
        return &m_cache[m_width * index(it)]; }
      const cache_type& operator()(const value_type& it, const unsigned entry ) const { 
        if( m_width * index(it) + entry > m_cache.size() ) ERROR("Invalid cache element: " << m_width * index(it) + entry > m_cache.size()  ); 
        return m_cache[m_width * index(it) + entry] ; }
      unsigned width()                                                          const { return m_width ; }
      
      template <typename functor_type> void set(const functor_type& functor, 
                                                unsigned int column,
                                                const std::string& key = "")
      { 
        for(const auto& element : *m_container) m_cache[ element.index() * m_width + column] = functor(element); 
        if( key != "" ) m_keys[column] = key; 
      }
      cache_type& operator()(const value_type& it, const unsigned entry ) { 
        auto pos = m_width * index(it) + entry;
        if( pos >= m_cache.size() ) ERROR("Out-of-bounds access: " << pos << " " << index(it) + entry);
        return m_cache[pos] ; }
      void setKey(const unsigned& column, const std::string& key ) { m_keys[column] = key ; }
      void print()
      {
        INFO( "width = " << m_width << ", size = " << m_size << " keys = " << vectorToString( m_keys , " ") << " cache size = " << m_cache.size() );
        for( unsigned int i = 0 ; i != m_width ; ++i ) std::cout << m_cache[i] << " ";
      }
    private:
      const container_type*     m_container; 
      std::vector<cache_type>   m_cache; 
      unsigned                  m_width; 
      unsigned                  m_size; 
      std::vector<std::string>  m_keys;
  };

}

#endif
