#include <array>

namespace AmpGen {
  template <class stored_type, class store_type, unsigned simd_size, bool modifiable> class scatter_iterator 
  {
    store_type*                        m_store; 
    std::array<stored_type, simd_size> m_buffer;
    size_t                   m_pos{0};
    public:
    int pos() const { return m_pos ; }
    scatter_iterator( const size_t& pos, store_type* store ) : 
      m_store(store), 
      m_pos(pos) { 
        if( m_store != nullptr && pos < m_store->aligned_size()) m_buffer = m_store->scatter(pos / simd_size );
      }
    stored_type* operator->() const { return &( m_buffer )[m_pos % simd_size]; }
    stored_type   operator*() const { return    m_buffer  [m_pos % simd_size]; }
    stored_type&  operator*()       { return    m_buffer  [m_pos % simd_size]; }
    scatter_iterator& operator++()  
    {
      m_pos++;
      if ( m_pos % simd_size == 0 )
      { 
        if constexpr(modifiable == true ) m_store->gather(m_buffer, (m_pos-1) / simd_size); 
        m_buffer = m_store->scatter( m_pos );
      }
      return *this;  
    }
    ~scatter_iterator()
    {
      if constexpr(modifiable == true)
      {
        if(m_store != nullptr && m_pos % simd_size != 0 ){
          m_store->gather(m_buffer, m_pos/simd_size); 
        }
      }
    }
    bool operator==( const scatter_iterator& rhs ) const { return m_pos == rhs.m_pos ; }
    bool operator!=( const scatter_iterator& rhs ) const { return m_pos != rhs.m_pos ; }
    friend int  operator-( const scatter_iterator& lhs, const scatter_iterator& rhs) { return lhs.pos() - rhs.pos() ; } 
  };
  template<unsigned simd_size, bool modifiable = false, class store_type> 
    auto make_scatter_iterator( const unsigned& pos, store_type* store) { 
      return scatter_iterator<typename store_type::value_type, store_type, simd_size, modifiable>(pos, store) ; }
}
