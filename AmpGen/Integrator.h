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

  template <size_t NROLL = 10>
    class Integrator
    {
      private:
        typedef const complex_t& arg;
        size_t                           m_counter = {0};
        std::array<Integral<arg>, NROLL> m_integrals;
        EventList*                       m_events  = {nullptr};
        void calculate()
        {
          integrateBlock();
          m_counter = 0;
        }
        void integrateBlock()
        {
          real_t re[NROLL] = {0};
          real_t im[NROLL] = {0};
          #pragma omp parallel for reduction(+: re, im)
          for ( size_t i = 0; i < m_events->size(); ++i ) {
            auto& evt = ( *m_events )[i];
            real_t w  = evt.weight() / evt.genPdf();
            for ( size_t roll = 0; roll < NROLL; ++roll ) {
              auto c = evt.getCache(m_integrals[roll].i) * std::conj(evt.getCache(m_integrals[roll].j));
              re[roll] += w * std::real(c);
              im[roll] += w * std::imag(c);
            }
          }
          real_t nv = m_events->norm();
          for ( size_t j = 0; j < m_counter; ++j )
            m_integrals[j].transfer( complex_t( re[j], im[j] ) / nv );
        }

      public:
        Integrator( EventList* events = nullptr ) : m_events( events ){}
        
        double sampleNorm()             { return m_events->norm(); } 
        bool isReady()            const { return m_events != nullptr; }
        EventList& events()             { return *m_events; } 
        const EventList& events() const { return *m_events; } 
        template <class T1, class T2>
        void addIntegral( const T1& f1, const T2& f2, const Integral<arg>::TransferFCN& tf )
        {
          addIntegralKeyed( m_events->getCacheIndex(f1), m_events->getCacheIndex(f2), tf );
        }
        void queueIntegral(const size_t& i, const size_t& j, complex_t* result){
          addIntegralKeyed(i, j, [result](arg& val){ *result = val ; } ); 
        }
        void queueIntegral(const size_t& c1, 
                           const size_t& c2, 
                           const size_t& i, 
                           const size_t& j, 
                           Bilinears* out, 
                           const bool& sim = true )
        {
          if( ! out->workToDo(i,j) )return;
          if( sim ) 
            addIntegralKeyed( c1, c2, [out,i,j]( arg& val ){ 
              out->set(i,j,val);
              if( i != j ) out->set(j,i, std::conj(val) ); } );
          else 
            addIntegralKeyed( c1, c2, [out,i,j]( arg& val ){ out->set(i,j,val); } );
        }
        void addIntegralKeyed(const size_t& c1, const size_t& c2, const Integral<arg>::TransferFCN& tf )
        {
          m_integrals[m_counter++] = Integral<arg>(c1, c2, tf);
          if ( m_counter == NROLL ) calculate();
        }

        void flush()
        {
          if ( m_counter == 0 ) return;
          calculate();
        }
        template <class EXPRESSION>
          void prepareExpression( const EXPRESSION& expression, const size_t& size_of = 0 )
          {
            if( m_events == nullptr ) return; 
            auto index = m_events->registerExpression( expression , size_of );
            m_events->updateCache( expression, index );
          }
        size_t getCacheIndex(const CompiledExpressionBase& expression) const {
          return m_events->getCacheIndex(expression);
        }
    };

  template <size_t NBINS = 100, size_t NROLL = 10>
    class BinnedIntegrator
    {
      private:
        typedef const std::vector<std::complex<double>>& arg;
        typedef std::function<void( arg )> TransferFCN;

        size_t                           m_counter = {0};
        std::vector<unsigned int>        m_view    = {0};
        std::vector<size_t>              m_slice   = {0};
        std::array<Integral<arg>, NROLL> m_integrals;
        EventList*                       m_events  = {nullptr};
        void calculate()
        {
          integrateBlock();
          m_counter = 0;
        }
        void integrateBlock()
        {
          double re[( NBINS + 1 ) * NROLL] = {0};
          double im[( NBINS + 1 ) * NROLL] = {0};
          auto ij                          = [&]( const Event& evt, const unsigned int& i, const unsigned int& j ) {
            return evt.getCache( i ) * std::conj( evt.getCache( j ) );
          };
          if ( m_slice.size() == 0 ) {
            #pragma omp parallel for reduction( + : re, im )
            for ( unsigned int i = 0; i < m_events->size(); ++i ) {
              auto& evt    = ( *m_events )[i];
              size_t binNo = m_view[i];
              double w     = evt.weight() / evt.genPdf();
              for ( unsigned int roll = 0; roll < NROLL; ++roll ) {
                auto c = ij( evt, m_integrals[roll].i, m_integrals[roll].j );
                DEBUG( "pos = " << roll * NBINS + binNo << " val = " << w * c );
                re[roll * NBINS + binNo] += w * std::real( c );
                im[roll * NBINS + binNo] += w * std::imag( c );
              }
            }
          } else {
            #pragma omp parallel for reduction( + : re, im )
            for ( unsigned int i = 0; i < m_slice.size(); ++i ) {
              auto& evt    = ( *m_events )[m_slice[i]];
              size_t binNo = m_view[i];
              double w     = evt.weight() / evt.genPdf();
              for ( unsigned int roll = 0; roll < NROLL; ++roll ) {
                auto c = ij( evt, m_integrals[roll].i, m_integrals[roll].j );
                re[roll * NBINS + binNo] += w * std::real( c );
                im[roll * NBINS + binNo] += w * std::imag( c );
              }
            }
          }
          double nv = m_events->norm();
          for ( size_t thisIntegral = 0; thisIntegral < m_counter; ++thisIntegral ) {
            std::vector<std::complex<double>> tmpBins( NBINS );
            size_t offset = thisIntegral * NBINS;
            for ( size_t nBin = 0; nBin < NBINS; ++nBin )
              tmpBins[nBin]   = std::complex<double>( re[offset + nBin], im[offset + nBin] ) / nv;
            m_integrals[thisIntegral].transfer( tmpBins );
          }
        }
      public:
        BinnedIntegrator( EventList* events = nullptr ) : m_events( events ) {}
        void setView( const std::function<size_t( const Event& )>& binNumber )
        {
          if ( m_slice.size() == 0 ) {
            if ( m_view.size() != m_events->size() ) m_view.resize( m_events->size() );
            for ( unsigned int i = 0; i < m_events->size(); ++i ) {

              m_view[i] = binNumber( ( *m_events )[i] );

              if ( m_view[i] >= NBINS ) {
                m_view[i] = NBINS;
                WARNING( "Event " << m_slice[i] << " bin number = " << m_view[i] << " is out of range!" );
              }
            }
          } else {
            if ( m_view.size() != m_slice.size() ) m_view.resize( m_slice.size() );
            for ( unsigned int i = 0; i < m_slice.size(); ++i ) {
              m_view[i] = binNumber( ( *m_events )[m_slice[i]] );
              if ( m_view[i] >= NBINS ) {
                m_view[i] = NBINS;
                WARNING( "Event " << m_slice[i] << " bin number = " << m_view[i] << " is out of range!" );
              }
            }
          }
        }
        void setSlice( const std::function<bool( const Event& )>& sliceFunction )
        {
          if ( m_slice.size() != 0 ) m_slice.clear();
          for ( unsigned int i = 0; i < m_events->size(); ++i )
            if ( sliceFunction( ( *m_events )[i] ) ) m_slice.push_back( i );
        }

        template <class T1, class T2>
          void addIntegral( const T1& f1, const T2& f2, const TransferFCN& tFunc )
          {
            m_integrals[m_counter++] = Integral<arg>( m_events->getCacheIndex( f1 ), m_events->getCacheIndex( f2 ), tFunc );
            if ( m_counter == NROLL ) calculate();
          }
        void flush()
        {
          if ( m_counter == 0 ) return;
          calculate();
        }
        template <class FCN>
          void update( FCN& fcn, std::array<Bilinears, NBINS>& normalisations )
          {
            auto mE   = fcn.matrixElements();
            auto size = mE.size();
            std::vector<size_t> toUpdate;
            std::vector<size_t> integralHasChanged( size * size );
            for ( size_t x = 0; x < size; ++x ) {
              auto& pdf = mE[x].amp;
              pdf.prepare();
              if ( !pdf.hasExternalsChanged() ) continue;
              m_events->updateCache( pdf, m_events->getCacheIndex( pdf ) );
              toUpdate.push_back( x );
            }
            for ( auto& i : toUpdate ) {
              DEBUG( "Updating: " << mE[i].decayTree->uniqueString() );
              for ( unsigned int j = 0; j < size; ++j ) {
                if ( integralHasChanged[i * size + j] ) continue;
                integralHasChanged[i * size + j] = true;
                integralHasChanged[j * size + i] = true;

                addIntegral( mE[i].amp, mE[j].amp,
                    [i, j, &normalisations]( const auto& val ) {
                    for ( unsigned int bin = 0; bin < NBINS; ++bin ) {
                    normalisations[bin].set( i, j, val[bin] );
                    if ( i != j ) normalisations[bin].set( j, i, std::conj( val[bin] ) );
                    }
                    } );
              }
            }
          }
    };
} // namespace AmpGen
#endif
