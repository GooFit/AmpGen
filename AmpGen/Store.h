#ifndef AMPGEN_STORE_H
#define AMPGEN_STORE_H

#include "AmpGen/simd/utils.h"
#include "AmpGen/EventList.h"
#include "AmpGen/ProfileClock.h"
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
      virtual ~Store() = default; 
      Store( const size_t& nEntries=0, const size_t& nFields=0) : 
        m_nEntries(nEntries), 
        m_nBlocks(utils::aligned_size<stored_type>( nEntries ) / utils::size<stored_type>::value ),
        m_nFields(nFields),
        m_store(m_nBlocks * m_nFields) {
          DEBUG("Calling basic constructor: " << m_nEntries << " nBlocks = " << m_nBlocks << " " << m_nFields << " " << nEntries ); 
        }

      template <typename functor_type> 
        void addFunctor( const functor_type& functor )
        {
          if( m_index.count(functor.name()) != 0 ) return; 
          unsigned vsize = functor.returnTypeSize() / sizeof(stored_type);
          DEBUG("Registering: " << functor.name() << " field: " << m_nFields << " " << functor.returnTypeSize() << " / " << sizeof(stored_type)  ); 
          std::vector<unsigned> offsets ( vsize ) ;
          for ( unsigned i = 0 ; i != vsize; ++i ) offsets[i] = m_nFields + i; 
          m_index[ functor.name() ] = offsets; 
          m_nFields += vsize;
        }

      template <typename functor_type> Store( const size_t& nEntries, const std::vector<functor_type>& functors)
      {
        allocate(nEntries, functors);
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type>
        Store( const size_t& nEntries, const functor_type& functor)
        {
          allocate( nEntries, {functor});
        }

      template <typename functor_type> void allocate( const size_t& nEntries, const std::vector<functor_type>& functors)
      {
        for(const auto& functor : functors) addFunctor( functor );
        resize( nEntries ); 
      }

      inline stored_type operator[]( const size_t& index ) const { return m_store[index]; }
      inline stored_type& operator[]( const size_t& index ) { return m_store[index]; }
      template <typename T> auto find( const T& t ) const { return m_index.find( t )->second; }

      inline size_t size()         const { return m_nEntries; }
      inline size_t size_raw()     const { return m_store.size(); }
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
      void store( const size_t& event0, const unsigned* index, const stored_type* item, const unsigned N = 1 )
      {
        for( unsigned i = 0 ; i != N; ++i ) (*this)(event0, index[i] ) = item[i]; 
      }

      template <typename functor_type> void update( const EventList& events, const functor_type& fcn )
      {  
      }
      void resize( std::size_t entries ) 
      {
        m_nEntries = entries;
        m_nBlocks  = utils::aligned_size<stored_type>(entries)/utils::size<stored_type>::value;
        DEBUG("resizing as: " << m_nEntries << " " << m_nBlocks << " " << m_nFields ); 
        m_store.resize(m_nBlocks * m_nFields); 
      }
    public:
      size_t m_nEntries{0};                                 ///< Number of entries, i.e. number of events 
      size_t m_nBlocks {0};                                 ///< Number of blocks, i.e. number of entries aligned to the size, divided by block size. 
      size_t m_nFields {0};                                 ///< Number of fields per entry 
      std::vector<stored_type> m_store;                     ///< Actual store of values  
      std::map<std::string, std::vector<unsigned>> m_index; ///< Index between name of functors and the stored values 
  };

  template <typename input_type, typename stored_type, Alignment align = SoA> class FunctionCache : public Store<stored_type, align> {
    public:
      FunctionCache() = default; 
      FunctionCache(const input_type* input, const std::size_t& nFields = 0 ) : 
        Store<stored_type, align>(m_input->size(), nFields ), m_input(input) {} 

      template <typename functor_type> FunctionCache(const input_type* input, const std::vector<functor_type>& functors)
      {
        allocate(input, functors); 
      }
      template <typename functor_type> void allocate( const input_type* input, const std::vector<functor_type>& functors ) 
      {
        m_input = input; 
        DEBUG( "Allocating: " << input->size() << " x " << functors.size() << " storage");
        for( const auto& f : functors ) Store<stored_type, align>::addFunctor(f);
        Store<stored_type, align>::resize( input->size() );  
        DEBUG("Finished allocation");
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type >
        void allocate( const input_type* input, const functor_type& functor ) 
        {
          m_input = input; 
          Store<stored_type, align>::addFunctor(functor);
          Store<stored_type, align>::resize( input->size() );  
          DEBUG("Finished allocation");
        }

      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type > FunctionCache( const input_type* input, const functor_type& functors)
        : FunctionCache( input, 1 )
      {
        Store<stored_type, align>::addFunctor(functors);
        Store<stored_type, align>::resize( input->size() );  
      }
      template <typename functor_type> void update(const functor_type& fcn)
      {
        const auto& f            = Store<stored_type, align>::find( fcn.name() ); 
        DEBUG("Updating: " << fcn.name() << " -> " << f[0] ); 
        if constexpr( ! std::is_same_v<input_type,EventList> )
        {
          auto stagger      = align == Alignment::AoS ? 1 : Store<stored_type,align>::m_nBlocks;
          auto fieldStagger = align == Alignment::AoS ? Store<stored_type, align>::m_nFields : 1;
          std::vector<size_t> offsets( f.size() );
          for( unsigned int i = 0 ; i != offsets.size(); ++i ) offsets[i] = f[i] * stagger; 
          if constexpr( std::is_same< typename functor_type::return_type, void >::value )
          {
            DEBUG("Evaluating bulk functor on: " << stagger << " " << m_input->nFields() << " " << fieldStagger  << " " << m_input->data()[0]  ); 
            fcn.batch(Store<stored_type, align>::aligned_size(), 
                m_input->nFields(), 
                fieldStagger,
                nullptr, 
                Store<stored_type, align>::data(), 
                offsets.data(),
                fcn.externBuffer().data(), 
                m_input->data());
            DEBUG( "Returning: " << real_v(Store<stored_type, align>::operator()(0,0).real()) ); //  << std::endl;
          }
          else 
          {
            DEBUG("Evaluating bulk function on [no-rto] : stagger:" << stagger << " Input Fields: " << m_input->nFields() << " field stagger: " << fieldStagger  << " input: " << m_input->data()[0] ); 

            fcn.batch(Store<stored_type,align>::aligned_size(), 
                m_input->nFields(), 
                fieldStagger, 
                Store<stored_type, align>::data() + f[0] *stagger, 
                fcn.externBuffer().data(), 
                m_input->data() ); 

          }
        }
        else { 
          auto p0 = f[0]; 
          auto s  = f.size(); 
          std::vector<size_t> offsets( s );
          std::iota( offsets.begin(), offsets.end(), 0 );
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t evt = 0; evt < m_input->size(); ++evt )
          { 
            if constexpr( std::is_same< typename functor_type::return_type, void >::value ) 
            {          
              std::vector<stored_type> buffer(s);
              fcn(buffer.data(), offsets.data(), fcn.externBuffer().data(), m_input->at(evt).address() ); 
              store(evt, f->second.data(), buffer.data(), s ); 
            } else {
              auto tmp = fcn( m_input->at(evt).address() );
              store( evt, f->second.data(), &tmp, s);
            }
          }
        }
      }
      const input_type*               m_input{nullptr};
  }; 
}
#if DEBUG_LEVEL ==1 
using aos_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::AoS>;
using soa_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::SoA>;

  ENABLE_DEBUG(aos_store)
ENABLE_DEBUG(soa_store)
#endif

#endif
