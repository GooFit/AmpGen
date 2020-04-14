#ifndef AMPGEN_INTEGRATOR_H
#define AMPGEN_INTEGRATOR_H

#include "AmpGen/Types.h"
#include "AmpGen/EventList.h"
#include <array>
#include <complex>
#include "AmpGen/simd/utils.h"
#include "AmpGen/Store.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventList.h"

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
    typedef std::function<void(const TYPE&)> TransferFCN;
    size_t i = {0};
    size_t j = {0};
    TransferFCN transfer;
    Integral() = default; 
    Integral(const size_t& i, const size_t& j, TransferFCN t) 
      : i(i), j(j), transfer(t) {}
  };
  
  class Integrator
  {
    typedef std::function<void(complex_t)> TransferFCN;

    public:
    Integrator() = default; 

    template <typename EventList_type, typename T> Integrator( const EventList_type* events, const std::vector<T>& expressions ={}, const size_t& size_of =0) : m_events(events)
    {
      if( events == nullptr ) {
        WARNING("No events specified, returning");
        return; 
      }
      m_cache = Store<complex_v, Alignment::SoA>(events->size(), expressions, size_of );
      m_weight.resize( events->nBlocks() );
      float_v norm_acc = 0.;
      for( size_t i = 0 ; i < events->nBlocks(); ++i )
      {
        m_weight[i]  = events->weight(i) / events->genPDF(i);
        norm_acc = norm_acc + m_weight[i]; 
      }
      m_norm = utils::sum_elements(norm_acc);
    }

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
    
    template <class return_type> return_type get( const unsigned& index, const unsigned& evt ) const ;
    template <class T> unsigned getCacheIndex( const T& t ) const { return m_cache.find(t) ; }
    double norm() const { return m_norm; }

    template <class T> void updateCache(const T& expression)
    {
      #if ENABLE_AVX2
      if( m_events != nullptr ) m_cache.update( static_cast<const EventListSIMD*>(m_events)->store(), expression );
      #else
      if( m_events != nullptr ) m_cache.update( static_cast<const EventList*>(m_events)->store(), expression );
      #endif
    }
    template <class T>
    const T* events() const { return static_cast<const T*>(m_events) ; }

    private:
    static constexpr size_t             N         = {8}; ///unroll factor
    size_t                              m_counter = {0};  ///
    std::array<Integral<complex_t>, N>  m_integrals;
    const void*                         m_events  = {nullptr};
    std::vector<float_v>                m_weight; 
    Store<complex_v, Alignment::SoA>    m_cache;      
    double                              m_norm    = {0};
    void integrateBlock();
  };

} // namespace AmpGen
#endif
