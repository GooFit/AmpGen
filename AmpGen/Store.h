#ifndef AMPGEN_STORE_H
#define AMPGEN_STORE_H

#include "AmpGen/simd/utils.h"
#include "AmpGen/EventList.h"
#ifdef _OPENMP 
#include <omp.h>
#endif

namespace AmpGen {

  enum Alignment {
    SoA, AoS
  };

  template <typename stored_type, Alignment align = SoA> class Store 
  {
    public:
      Store( const size_t& nEntries=0, const size_t& nFields=0) : 
        m_nEntries(nEntries), 
        m_nBlocks(utils::aligned_size<stored_type>( nEntries ) / utils::size<stored_type>::value ),
        m_nFields(nFields),
        m_store(m_nBlocks * m_nFields) {}
      
      template <typename functor_type> 
      void addFunctor( const functor_type& functor, unsigned fieldsPerFunctor=0 )
      {
        if( m_index.count(functor.name()) == 0 )
        {
          auto vsize = fieldsPerFunctor == 0 ? functor.returnTypeSize() / sizeof(stored_type) : fieldsPerFunctor;
          DEBUG("Registering: " << functor.name() << " field = " << m_nFields ); 
          m_index[ functor.name() ] = std::make_pair(m_nFields, vsize);
          m_nFields += vsize;
        }
      }
      template <typename functor_type> void allocate( const size_t& nEntries, const std::vector<functor_type>& functors, const size_t& fieldsPerFunctor = 0)
      {
        for(const auto& functor : functors) addFunctor( functor, fieldsPerFunctor);
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_store.resize(m_nBlocks * m_nFields); 
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type > 
        void allocate( const size_t& nEntries, const functor_type& functor, const size_t& fieldsPerFunctor = 0)
      {
        addFunctor(functor, fieldsPerFunctor);
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_store.resize(m_nBlocks * m_nFields); 
      }
       
      template <typename functor_type> Store( const size_t& nEntries, const std::vector<functor_type>& functors, const size_t& fieldsPerFunctor = 0)
      {
        allocate(nEntries, functors, fieldsPerFunctor);
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type>
      Store( const size_t& nEntries, const functor_type& functor, const size_t& fieldsPerFunctor=0 )
      {
        allocate( nEntries, {functor}, fieldsPerFunctor);
      }

      inline stored_type operator[]( const size_t& index ) const { return m_store[index]; }
      inline stored_type& operator[]( const size_t& index ) { return m_store[index]; }
      template <typename T> unsigned find( const T& t ) const { return m_index.find( t.name() )->second.first; }

      inline size_t size()         const { return m_nEntries; }
      inline size_t nBlocks()      const { return m_nBlocks; }
      inline size_t nFields()      const { return m_nFields; }
      inline size_t aligned_size() const { return m_nBlocks * utils::size<stored_type>::value ; }
      inline const stored_type& operator()(const size_t& index, const size_t& field) const
      {
        if constexpr( align == Alignment::SoA ) return m_store[ field * m_nBlocks + index] ; 
        else return m_store[index*m_nFields+field]; 
      }
      template <typename return_type> 
      inline const return_type get(const size_t& index, const size_t& field ) const 
      {
        return utils::at( operator()( index / utils::size<stored_type>::value, field ), index % utils::size<stored_type>::value );
      }
      inline const stored_type* data() const { return m_store.data(); }
      inline stored_type* data() { return m_store.data() ;}
      inline stored_type& operator()(const size_t& index, const size_t& field)
      {
        if constexpr( align == Alignment::SoA ) return m_store[ field * m_nBlocks + index] ; 
        else return m_store[index*m_nFields+field]; 
      }
      
      void resize(const size_t& nEntries, const size_t& nFields )
      {
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_nFields  = nFields;
        m_store.resize(m_nBlocks * m_nFields); 
        m_index.clear();
      }
      void clear() { m_store.clear(); m_index.clear() ; }
      void store( const size_t& event0, const size_t& index0, const stored_type* item, const unsigned N = 1 )
      {
        if constexpr( align == Alignment::AoS ) 
          std::memcpy( &(*this)(event0, index0) , item, N * sizeof( stored_type ) );
        else 
        {
          for( unsigned i = 0 ; i != N ; ++i ) (*this)(event0, index0 +i ) = item[i];
        }   
      }

      template <typename functor_type, typename input_type> void update(const Store<input_type, Alignment::AoS>& is, const functor_type& fcn)
      {
        auto f = m_index.find( fcn.name() ); 
        if( f == m_index.end() ) FATAL("Expression: " << fcn.name() << " is not registed");
        auto [p0, s] = f->second;
        DEBUG("Updating: " << fcn.name() << " index = " << p0 << " size_of = " << s << " on store: " << is.size() << " blocks = " << is.nBlocks() << " fields = " << is.nFields () ); 
        if constexpr( align == Alignment::AoS )
        {
          if constexpr( std::is_same< typename functor_type::return_type, void >::value )
          fcn.batch(aligned_size(), is.nFields(), m_nFields, nullptr,  m_store.data() + p0, 1, fcn.externBuffer().data(), is.data());  
          if constexpr( ! std::is_same< typename functor_type::return_type, void >::value )
          fcn.batch(aligned_size(), is.nFields(), m_nFields         , m_store.data()  + p0   , fcn.externBuffer().data(), is.data());
        }
        else 
        {
          if constexpr( std::is_same< typename functor_type::return_type, void >::value)
            fcn.batch(aligned_size(), is.nFields(), 1, nullptr, m_store.data() + p0*m_nBlocks, m_nBlocks, fcn.externBuffer().data(), is.data() ); 
          else 
            fcn.batch(aligned_size(), is.nFields(), 1         , m_store.data() + p0*m_nBlocks           , fcn.externBuffer().data(), is.data() ); 
        }
      }
      template <typename functor_type> void update( const EventList& events, const functor_type& fcn )
      {  
        auto f = m_index.find( fcn.name() ); 
        if( f == m_index.end() ) FATAL("Expression: " << fcn.name() << " is not registed");
        //auto& [p0, s] = f->second; /// bug in the C++ standard. Such fun. 
        auto p0 = f->second.first;
        auto s  = f->second.second; 
        if constexpr( std::is_same< typename functor_type::return_type, void >::value ) 
        {
          
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t evt = 0; evt < events.size(); ++evt )
          { 
            std::vector<stored_type> buffer(s);
            fcn(buffer.data(), 1, fcn.externBuffer().data(), events[evt].address() ); 
            store(evt, p0, buffer.data(), s ); 
          }
        }
        else {
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t evt = 0; evt < events.size(); ++evt )
          {
            auto tmp = fcn( events[evt].address() );
            store( evt, p0, &tmp, s);
          }
        }
      }

    private:
      size_t m_nEntries{0}; /// Number of entries, i.e. number of events 
      size_t m_nBlocks {0}; /// Number of blocks, i.e. number of entries aligned to the size, divided by block size. 
      size_t m_nFields {0}; /// Number of fields per entry 
      std::vector<stored_type> m_store; 
      std::map<std::string, std::pair<unsigned, unsigned>> m_index; 
  };
}
// using aos_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::AoS>;
// using soa_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::SoA>;
// 
// ENABLE_DEBUG(aos_store)
// ENABLE_DEBUG(soa_store)
#endif
