#ifndef AMPGEN_INTEGRATOR2_H
#define AMPGEN_INTEGRATOR2_H 1

#include "AmpGen/Integrator.h"

namespace AmpGen {

  class Integrator2
  {
      typedef const complex_t& arg;
      typedef std::function<void(arg)> TransferFCN;

    public:
      explicit Integrator2( const EventList* events = nullptr );

      bool isReady()            const;
      const EventList& events() const;
      void queueIntegral(const size_t& c1, 
          const size_t& c2, 
          const size_t& i, 
          const size_t& j, 
          Bilinears* out, 
          const bool& sim = true);
      void addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc );
      void queueIntegral(const size_t& i, const size_t& j, complex_t* result);
      void flush();
      void setBuffer( complex_t* pos, const complex_t& value, const size_t& size );
      void setBuffer( complex_t* pos, const std::vector<complex_t>& value, const size_t& size);
      
      template <class T> size_t getCacheIndex(const T& expression) const 
      {
        return m_index.find(expression.name())->second;
      }
      template <class T> void prepareExpression( const T& expression, const size_t& size_of = 0 )
      {
        if( m_events == nullptr ) return; 
        size_t vsize = size_of == 0 ? expression.returnTypeSize() / sizeof(complex_t) : size_of;
        auto it = m_index.find( expression.name() );
        auto index = 0;
        if( it == m_index.end() )
        {
          index = m_buffer.size();
          m_index[ expression.name() ] = index;
          m_buffer.resize(index+vsize);
          for(size_t j = 0 ; j != vsize; ++j )
            m_buffer[index+j].resize( m_events->size() );
        }
        else index = it->second; 
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for ( size_t i = 0; i < m_events->size(); ++i )
        {
          auto v = expression(m_events->at(i).address());
          setBuffer( &(m_buffer[index][i]), v, vsize );
        }
      }
    
    private:
      static constexpr size_t             N         = {10}; ///unroll factor
      size_t                              m_counter = {0};  ///
      std::array<Integral<arg>, N>        m_integrals;
      const EventList*                    m_events  = {nullptr};
      std::vector<std::vector<complex_t>> m_buffer;
      std::vector<double>                 m_weight; 
      std::map<std::string, size_t>       m_index; 
      double                              m_norm    = {0};
      void integrateBlock();
  };
}
#endif
