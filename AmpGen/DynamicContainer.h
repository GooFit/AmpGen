#ifndef AMPGEN_DYNAMICCONTAINER_H
#define AMPGEN_DYNAMICCONTAINER_H
#include "AmpGen/MsgService.h"
#include <functional>
#include <vector>

namespace AmpGen
{
  /**@class DynamicContainer
   *  Template class for a DynamicContainer, containing a temporary buffer of objects and
   *  a generator function that refreshes the buffer once the current buffer is exhausted.
   *  Can only be navigated using quasi foward-iterators.
   *  Designed such that very large sets can be operated on if only some value from the set is 
   *  required at the end of the computation, for example, 
   *  for computing large sets of multidimensional integral is is more efficient to tranch 
   *  the computation onto blocks of events.
   */
  template <class TYPE, class CONTAINER_TYPE = std::vector<TYPE>>
  class DynamicContainer
  {
  public:
    class Iterator
    {
    private:
      DynamicContainer* m_container;
      size_t m_pos;
    public:
      Iterator( const size_t& pos, DynamicContainer* parent ) : m_container( parent ), m_pos( pos ) {}
      TYPE* operator->() const { return &( ( *m_container )[m_pos] ); }
      TYPE operator*() const { return ( *m_container )[m_pos]; }
      Iterator operator++() const
      {
        if ( m_pos + 1 % m_container->m_buffer.size() == 0 && m_pos + 1 != m_container->m_totalSize ) {
          m_container->refresh();
        }
        return Iterator( m_pos + 1, m_container );
      }
      Iterator& operator++()
      {
        m_pos++;
        if ( m_pos % m_container->m_buffer.size() == 0 && m_pos % m_container->m_totalSize ) {
          m_container->refresh();
        }
        return *this;
      }
      bool operator==( const Iterator& rhs ) const { return m_pos == rhs.m_pos; }
      bool operator!=( const Iterator& rhs ) const { return m_pos != rhs.m_pos; }
    };
  private:
    CONTAINER_TYPE m_buffer;
    size_t m_totalSize;
    std::function<void( CONTAINER_TYPE& )> m_generator;
  public:
    template <class GENERATOR> DynamicContainer( const size_t& totalSize, const GENERATOR& generator ) : 
      m_totalSize( totalSize ), m_generator( generator )
      {
        refresh();
      }
    Iterator begin() { return Iterator( 0, this ); }
    Iterator end() { return Iterator( m_totalSize, this ); }
    const size_t& size() const { return m_totalSize; }
    TYPE& operator[]( const size_t& pos ) { return m_buffer[pos % m_buffer.size()]; }
    const TYPE& operator[]( const size_t& pos ) const { return m_buffer[pos % m_buffer.size()]; }
    void refresh()
    {
      m_buffer.clear();
      m_generator( m_buffer );
    }
    CONTAINER_TYPE& buffer() { return m_buffer; }
  };
} // namespace AmpGen
#endif
