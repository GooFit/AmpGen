#ifndef AMPGEN_INTEGRATOR_H
#define AMPGEN_INTEGRATOR_H

#include "AmpGen/Types.h"
#include "AmpGen/EventList.h"
#include "AmpGen/CompiledExpressionBase.h"
#include <array>
#include <complex>

/*
 *  Calculates Bilinears A_i A_j^* integrated over the phase-space.
 *  Integrates in blocks of (i,j) such that integrals can be queued and evaluated in blocks
 *  to optimise cache throughput.
 */

namespace AmpGen
{
  class Bilinears 
  {
    private: 
      size_t rows;
      size_t cols;
      std::vector<complex_t> norms;
      std::vector<bool> markAsZero; 
      std::vector<bool> calculate;   
    public:
      Bilinears( const size_t& r = 0, const size_t& c = 0 );
      complex_t get(const size_t& x, const size_t& y) const;
      template <class T>
      complex_t get(const size_t& x, const size_t& y, T* integ = nullptr, const size_t& kx=0, const size_t& ky=0){
        if( integ != nullptr ) integ->queueIntegral(kx, ky, &norms[x*cols+y]);
        /// will return the wrong answer for now, but queues for later.. 
        return norms[x*cols+y];
      }
      void      set(const size_t& x, const size_t& y,  const complex_t& f );
      void  setZero(const size_t& x, const size_t& y);
      void resetCalculateFlags();
      complex_t& operator()( const size_t& x, const size_t& y );
      bool   isZero(const size_t& x, const size_t& y);
      bool workToDo(const size_t& x, const size_t& y) const;
      void   resize(const size_t& r, const size_t& c = 1 );
  };

  template <class TYPE = complex_t> struct Integral 
  {
    typedef std::function<void(TYPE)> TransferFCN;
    size_t i = {0};
    size_t j = {0};
    TransferFCN transfer;
    Integral() = default; 
    Integral(const size_t& i, const size_t& j, TransferFCN t) 
      : i(i), j(j), transfer(t) {}
  };
  class Integrator
  {
      typedef const complex_t& arg;
      typedef std::function<void(arg)> TransferFCN;

    public:
      explicit Integrator( const EventList* events = nullptr );

      bool isReady()            const;
      const EventList* events() const;
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
      complex_t get(const unsigned& i, const unsigned& evt) const { return m_cache[i * m_events->size() + evt ]; }
      template <class T> unsigned getCacheIndex( const T& t ) const { return m_index.find( t.name() )->second.first; }
      double norm() const { return m_norm; }
      template <class T> void allocate( const std::vector<T>& expressions, const size_t& size_of = 0)
      {
        if( m_events == nullptr ) return; 
        unsigned totalSize = 0; 
        for( unsigned i = 0; i != expressions.size(); ++i ){
          size_t vsize = size_of == 0 ? expressions[i].returnTypeSize() / sizeof(complex_t) : size_of;
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
        INFO("Preparing: " << expression.name() << " index = " << p0 << " with: " << s << " values" );
        if constexpr( std::is_same< typename T::return_type, void >::value ) 
        {
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t i = 0; i < m_events->size(); ++i )
          { 
            std::vector<complex_t> buf(s);
            expression(&buf[0], expression.externBuffer().data(), m_events->at(i).address() );
            for( unsigned j = 0; j != s; ++j ) m_cache[ (p0+j) * m_events->size() + i] = buf[j];
          }
        }
        else {
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t i = 0; i < m_events->size(); ++i )
            setBuffer( &(m_cache[p0 * m_events->size() +i] ), expression(m_events->at(i).address()),s );
        }
      }
    
    private:
      static constexpr size_t             N         = {10}; ///unroll factor
      size_t                              m_counter = {0};  ///
      std::array<Integral<arg>, N>        m_integrals;
      const EventList*                    m_events  = {nullptr};
      std::vector<complex_t>              m_cache;
      std::vector<double>                 m_weight; 
      std::map<std::string, std::pair<unsigned, unsigned>>       m_index; 
      double                              m_norm    = {0};
      void integrateBlock();
  };
} // namespace AmpGen
#endif
