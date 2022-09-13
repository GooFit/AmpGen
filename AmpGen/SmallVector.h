#include <array>

namespace AmpGen {
  template <typename type, unsigned max_size> 
    struct SmallVector {
      typedef type value_type; 

      std::array<type, max_size> store = {0}; 
      unsigned size = 0; 
      auto begin() const { return store.begin(); }
      auto end()   const { return store.begin() + size; } 
      auto begin()  { return store.begin(); }
      auto end()    { return store.begin() + size; } 
      value_type& operator[]       (unsigned i)       { return store[i]; } 
      const value_type& operator[] (unsigned i) const { return store[i]; } 
      void push_back( const value_type& thing ){ store[size++] = thing; }
      template <typename iterator_type, typename other_iterator_type> 
        void insert(iterator_type mbegin, other_iterator_type ibegin, other_iterator_type iend )
        {  
          for( auto it = ibegin; it != iend; ++it ){
            *( mbegin + (it - ibegin) ) = *it;  
          }    
          size += (iend - ibegin ); 
        }
      SmallVector() = default; 
      SmallVector( std::initializer_list<value_type>&& values )
      {
        insert( begin(), values.begin(), values.end() );
      }
      bool operator==( const SmallVector& other ) const
      {
        return other.store == this->store; 
      }
    };
}
