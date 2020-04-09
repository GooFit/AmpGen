#ifndef AMPGEN_INTEGRATORSIMD_H
#define AMPGEN_INTEGRATORSIMD_H 1

#if ENABLE_AVX2

#include "AmpGen/Integrator.h"
#include "AmpGen/simd/avx2_types.h"
#include "AmpGen/EventListSIMD.h"

namespace AmpGen {
  /// test /// 
  class IntegratorSIMD
  {
      typedef const complex_t& arg;
      typedef std::function<void(arg)> TransferFCN;

    public:
      explicit IntegratorSIMD( const EventListSIMD* events = nullptr );

      bool isReady()            const;
      void queueIntegral(const size_t& c1, 
          const size_t& c2, 
          const size_t& i, 
          const size_t& j, 
          Bilinears* out, 
          const bool& sim = true);
      void addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc );
      void queueIntegral(const size_t& i, const size_t& j, complex_t* result);
      void flush();
      complex_t get(const unsigned& i, const unsigned& evt) const { return m_cache[i * m_events->size() + evt/float_v::size ].at( evt % float_v::size ); }
      template <class T> unsigned getCacheIndex( const T& t ) const { return m_index.find( t.name() )->second.first; }
      double norm() const { return m_norm; }
      const EventListSIMD* events() const { return m_events; }
      template <class T> void allocate( const std::vector<T>& expressions, const size_t& size_of = 0)
      {
        if( m_events == nullptr ) return; 
        unsigned totalSize = 0; 
        for( unsigned i = 0; i != expressions.size(); ++i ){
          size_t vsize = size_of == 0 ? expressions[i].returnTypeSize() / sizeof(AVX2::complex_t) : size_of;
          m_index[ expressions[i].name() ] = std::make_pair(totalSize, vsize);
          totalSize += vsize;
        }
        m_cache.resize( m_events->size() * totalSize );
      }
      
      template <class T> void prepareExpression(const T& expression)
      {
        if( m_events == nullptr ) return; 
        auto f = m_index.find( expression.name() ); 
        if( f == m_index.end() ) FATAL("Expression: " << expression.name() << " is not registed");
        auto [p0, s] = f->second;
        expression.batch(m_cache.data() + p0*m_events->aligned_size(), 
                         m_events->aligned_size(), 
                         m_events->eventSize(), 
                         1, 
                         expression.externBuffer().data(), 
                         m_events->data() ); 
      }
    
    private:
      static constexpr size_t             N         = {10}; ///unroll factor
      size_t                              m_counter = {0};  ///
      std::array<Integral<arg>, N>        m_integrals;
      const EventListSIMD*                   m_events  = {nullptr};
      std::vector<complex_v>              m_cache;
      std::vector<float_v>                m_weight; 
      std::map<std::string, std::pair<unsigned, unsigned>>       m_index; 
      double                                                     m_norm    = {0};
      void integrateBlock();
  };
}
#endif
#endif
