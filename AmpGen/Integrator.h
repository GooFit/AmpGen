#ifndef AMPGEN_INTEGRATOR_H
#define AMPGEN_INTEGRATOR_H

#include "AmpGen/Types.h"
#include "AmpGen/EventList.h"
#include <array>
#include <complex>

/*
 *  Calculates Bilinears A_i A_j^* integrated over the phase-space.
 *  Integrates in blocks of (i,j) such that integrals can be queued and evaluated in blocks
 *  to optimise cache throughput.
 */

namespace AmpGen
{
  struct Bilinears {
    size_t rows;
    size_t cols;
    std::vector<complex_t> norms;
    Bilinears( const size_t r = 0, const size_t c = 0 ) : rows( r ), cols( c ), norms( c * r ) {}

    complex_t get( const unsigned int& x, const unsigned int& y ) const { return norms[x * cols + y]; }
    void set( const size_t& x, const size_t& y, const complex_t& f ) { norms[x * cols + y] = f; }
    void resize( const size_t r, const size_t c = 1 )
    {
      rows = r;
      cols = c;
      norms.resize( r * c );
    }
  };

  template <class TYPE = complex_t>
  struct Integral {
    typedef std::function<void( TYPE )> TransferFCN;
    unsigned int i;
    unsigned int j;
    TransferFCN transfer;
    Integral( const unsigned int& _i, const unsigned int& _j, TransferFCN _t ) : i( _i ), j( _j ), transfer( _t ) {}
    Integral() : i( 0 ), j( 0 ) {}
  };

  template <size_t NROLL = 10>
  class Integrator
  {

  public:
    EventList* sim;

  private:
    typedef const complex_t& arg;
    typedef std::function<void( arg )> TransferFCN;

    size_t m_counter;
    std::array<Integral<arg>, NROLL> m_integrals;
    void calculate()
    {
      integrateBlock();
      m_counter = 0;
    }
    void integrateBlock()
    {
      real_t re[NROLL] = {0};
      real_t im[NROLL] = {0};
      auto ij          = [&]( const Event& evt, const size_t& i, const size_t& j ) {
        return evt.getCache( i ) * std::conj( evt.getCache( j ) );
      };
      #pragma omp parallel for reduction( + : re, im )
      for ( size_t i = 0; i < sim->size(); ++i ) {
        auto& evt = ( *sim )[i];
        real_t w  = evt.weight() / evt.genPdf();
        for ( size_t roll = 0; roll < NROLL; ++roll ) {
          auto c = ij( evt, m_integrals[roll].i, m_integrals[roll].j );
          re[roll] += w * std::real( c );
          im[roll] += w * std::imag( c );
         // INFO( "Adding to roll " << roll << " = " << m_integrals[roll].i << " " << m_integrals[roll].j << " " << evt.getCache( m_integrals[roll].i ) << " x " << evt.getCache( m_integrals[roll].j ) );
        }
      }
      real_t nv = sim->norm();
      for ( size_t j = 0; j < m_counter; ++j )
        m_integrals[j].transfer( complex_t( re[j], im[j] ) / nv );
    }
  public:
    Integrator( EventList* _sim = nullptr ) : sim( _sim ), m_counter( 0 ) {}

    template <class T1, class T2>
    void addIntegral( const T1& f1, const T2& f2, const TransferFCN& tFunc )
    {
      addIntegralKeyed( sim->getCacheIndex(f1), sim->getCacheIndex(f2), tFunc );
    }
    void addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc )
    {
      m_integrals[m_counter++] = Integral<arg>(c1,c2,tFunc);
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
      auto index = sim->registerExpression( expression , size_of );
      sim->updateCache( expression, index );
    }
  };

  template <size_t NBINS = 100, size_t NROLL = 10>
  class BinnedIntegrator
  {
  public:
    EventList* sim;

  private:
    typedef const std::vector<std::complex<double>>& arg;
    typedef std::function<void( arg )> TransferFCN;

    size_t m_counter;
    std::vector<unsigned int> m_view;
    std::vector<size_t> m_slice;
    std::array<Integral<arg>, NROLL> m_integrals;
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
        for ( unsigned int i = 0; i < sim->size(); ++i ) {
          auto& evt    = ( *sim )[i];
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
          auto& evt    = ( *sim )[m_slice[i]];
          size_t binNo = m_view[i];
          double w     = evt.weight() / evt.genPdf();
          for ( unsigned int roll = 0; roll < NROLL; ++roll ) {
            auto c = ij( evt, m_integrals[roll].i, m_integrals[roll].j );
            re[roll * NBINS + binNo] += w * std::real( c );
            im[roll * NBINS + binNo] += w * std::imag( c );
          }
        }
      }
      double nv = sim->norm();
      for ( size_t thisIntegral = 0; thisIntegral < m_counter; ++thisIntegral ) {
        std::vector<std::complex<double>> tmpBins( NBINS );
        size_t offset = thisIntegral * NBINS;
        for ( size_t nBin = 0; nBin < NBINS; ++nBin )
          tmpBins[nBin]   = std::complex<double>( re[offset + nBin], im[offset + nBin] ) / nv;
        m_integrals[thisIntegral].transfer( tmpBins );
      }
    }

  public:
    void setView( const std::function<size_t( const Event& )>& binNumber )
    {
      if ( m_slice.size() == 0 ) {
        if ( m_view.size() != sim->size() ) m_view.resize( sim->size() );
        for ( unsigned int i = 0; i < sim->size(); ++i ) {

          m_view[i] = binNumber( ( *sim )[i] );

          if ( m_view[i] >= NBINS ) {
            m_view[i] = NBINS;
            WARNING( "Event " << m_slice[i] << " bin number = " << m_view[i] << " is out of range!" );
          }
        }
      } else {
        if ( m_view.size() != m_slice.size() ) m_view.resize( m_slice.size() );
        for ( unsigned int i = 0; i < m_slice.size(); ++i ) {
          m_view[i] = binNumber( ( *sim )[m_slice[i]] );
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
      for ( unsigned int i = 0; i < sim->size(); ++i )
        if ( sliceFunction( ( *sim )[i] ) ) m_slice.push_back( i );
    }
    BinnedIntegrator( EventList* _sim = nullptr ) : sim( _sim ), m_counter( 0 ), m_slice( 0 ) {}

    template <class T1, class T2>
    void addIntegral( const T1& f1, const T2& f2, const TransferFCN& tFunc )
    {
      m_integrals[m_counter++] = Integral<arg>( sim->getCacheIndex( f1 ), sim->getCacheIndex( f2 ), tFunc );
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
        auto& pdf = mE[x].pdf;
        pdf.prepare();
        if ( !pdf.hasExternalsChanged() ) continue;
        sim->updateCache( pdf, sim->getCacheIndex( pdf ) );
        toUpdate.push_back( x );
      }
      for ( auto& i : toUpdate ) {
        DEBUG( "Updating: " << mE[i].decayTree->uniqueString() );
        for ( unsigned int j = 0; j < size; ++j ) {
          if ( integralHasChanged[i * size + j] ) continue;
          integralHasChanged[i * size + j] = true;
          integralHasChanged[j * size + i] = true;

          addIntegral( mE[i].pdf, mE[j].pdf,

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
