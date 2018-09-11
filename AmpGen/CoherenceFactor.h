#ifndef AMPGEN_COHERENCEFACTOR_H
#define AMPGEN_COHERENCEFACTOR_H
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "AmpGen/BinDT.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

namespace AmpGen
{
  struct CoherenceEvent {
    std::array<double, 5> values;
    std::complex<double> amp1;
    std::complex<double> amp2;
    double weight;
    uint64_t flag;
  };

  struct CoherenceCalculator {
    FastCoherentSum* m_pdf1;
    FastCoherentSum* m_pdf2;
    Bilinears R;
    Bilinears N1;
    Bilinears N2;
    double rTransform;
    double phiTransform;
    CoherenceCalculator() = default;
    CoherenceCalculator( FastCoherentSum* pdf1, FastCoherentSum* pdf2 )
        : m_pdf1( pdf1 )
        , m_pdf2( pdf2 )
        , R( m_pdf1->size(), m_pdf2->size() )
        , N1( m_pdf1->size(), m_pdf1->size() )
        , N2( m_pdf2->size(), m_pdf2->size() )
        , rTransform( 1 )
        , phiTransform( 0 )
    {
    }

    void setTransform( const double& rT, const double& pT )
    {
      rTransform   = rT;
      phiTransform = pT;
    }
    double getN2() const
    {
      double accN2 = 0;
      for ( unsigned int i = 0; i < m_pdf2->size(); ++i ) {
        accN2 += std::norm( ( *m_pdf2 )[i].coefficient ) * std::abs( N2.get( i, i ) );
        for ( unsigned int j = i + 1; j < m_pdf2->size(); ++j ) {
          accN2 +=
              2 * std::real( ( *m_pdf2 )[i].coefficient * std::conj( ( *m_pdf2 )[j].coefficient ) * N2.get( i, j ) );
        }
      }
      return accN2;
    }

    double getN1() const
    {
      double accN1 = 0;
      for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
        accN1 += std::norm( ( *m_pdf1 )[i].coefficient ) * std::abs( N1.get( i, i ) );
        for ( unsigned int j = i + 1; j < m_pdf1->size(); ++j ) {
          accN1 +=
              2 * std::real( ( *m_pdf1 )[i].coefficient * std::conj( ( *m_pdf1 )[j].coefficient ) * N1.get( i, j ) );
        }
      }
      return accN1 * rTransform * rTransform;
    }

    std::complex<double> getVal()
    {
      std::complex<double> accR( 0, 0 );
      m_pdf1->transferParameters();
      m_pdf2->transferParameters();
      for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
        for ( unsigned int j = 0; j < m_pdf2->size(); ++j ) {
          accR += ( *m_pdf1 )[i].coefficient * std::conj( ( *m_pdf2 )[j].coefficient ) * R.get( i, j );
        }
      }
      return std::complex<double>( cos( phiTransform ), sin( phiTransform ) ) * rTransform * accR /
             sqrt( getN1() * getN2() );
    }
  };

  struct HadronicParameters {

    HadronicParameters() : n1( 0 ), n2( 0 ), coherence( 0, 0 ), nFills( 0 ), id( 0 ), binID( 0 ) {}

    double n1;
    double n2;
    std::complex<double> coherence;
    unsigned int nFills;
    unsigned int id;
    unsigned int binID;
    std::complex<double> getCoherence() const { return coherence / sqrt( n1 * n2 ); }
    void reset()
    {
      n1        = 0;
      n2        = 0;
      coherence = 0;
      nFills    = 0;
      id        = 0;
    }
    void scale_p1( const double& sf )
    {
      n1 /= sf;
      coherence /= sqrt( sf );
    }
    void scale_p2( const double& sf )
    {
      n2 /= sf;
      coherence /= sqrt( sf );
    }

    void rotate( const double& angle ) { coherence *= std::complex<double>( cos( angle ), sin( angle ) ); }
    double d() const { return -std::arg( coherence ); }
    double R() const { return std::abs( getCoherence() ); }
    double r() const { return sqrt( n1 / n2 ); }
    void clear() { n1 = ( 0 ), n2 = ( 0 ), coherence = std::complex<double>( 0, 0 ); }

    void add( const std::complex<double>& f1, const std::complex<double>& f2, const double& weight = 1 )
    {
      n1 += weight * std::norm( f1 );
      n2 += weight * std::norm( f2 );
      coherence += weight * f1 * std::conj(f2) ;
      nFills++;
    }
    HadronicParameters( const double& _R, const double& _d, const double& k1 )
        : n1( k1 * k1 ), n2( 1 ), coherence( k1 * _R * cos( M_PI * _d / 180. ), k1 * _R * sin( M_PI * _d / 180. ) ){};
    HadronicParameters( const double& _R, const double& _d, const double& k1, const double& k2 )
        : n1( k1 ), n2( k2 ), coherence( sqrt( k1 * k2 ) * _R * cos( _d ), sqrt( k1 * k2 ) * _R * sin( _d ) ){};

    HadronicParameters operator+=( const HadronicParameters& other )
    {
      n1 += other.n1;
      n2 += other.n2;
      coherence += other.coherence;
      nFills += other.nFills;
      return *this;
    }
    HadronicParameters operator-=( const HadronicParameters& other )
    {
      n1 -= other.n1;
      n2 -= other.n2;
      coherence -= other.coherence;
      nFills -= other.nFills;
      return *this;
    }

    std::string to_string() const
    {
      auto c = getCoherence();

      std::string total = "R = " + round( R() , 3 ) + " δ = " + round( 180 * d() / M_PI, 2 ) +
                          " r = " + round( r(), 4 ) + " K = " + round( n1, 6 ) + " K' = " + round( n2, 6 ) + " N = " +
                          std::to_string( nFills ) + " ci = " + std::to_string( std::real( c ) ) + " si = " +
                          std::to_string( std::imag( c ) );
      return total;
    }
  };
  HadronicParameters operator+( const HadronicParameters& a, const HadronicParameters& b )
  {
    HadronicParameters hp;
    hp.n1        = a.n1 + b.n1;
    hp.n2        = a.n2 + b.n2;
    hp.coherence = a.coherence + b.coherence;
    hp.nFills    = a.nFills + b.nFills;
    return hp;
  }

  class CoherenceFactor
  {
  private:
    HadronicParameters m_global;
    double m_globalPhase;
    double m_globalR;
    FastCoherentSum* m_pdf1;
    FastCoherentSum* m_pdf2;
    BinDT m_voxels;
    std::vector<HadronicParameters> m_paramsPerVoxel;
    std::map<unsigned int, unsigned int> m_nodeID2Bin;
    unsigned int m_nBins;

    CoherenceCalculator m_calc;
    EventList* m_data;
    EventType m_type;

  public:
    BinDT& getDT() { return m_voxels; }
    double phase();
    double getR();
    std::vector<size_t> getNumberOfEventsInEachBin( const EventList& events ) const;
    HadronicParameters calculateGlobalCoherence( EventList* data = nullptr );
    void writeToFile( const std::string& filename );

    HadronicParameters getVal() const;
    CoherenceFactor( const std::string& filename );
    CoherenceFactor();
    CoherenceFactor( FastCoherentSum* pdf1, FastCoherentSum* pdf2, EventList* data = nullptr );

    CoherenceFactor( FastCoherentSum* pdf1, FastCoherentSum* pdf2, const EventType& type );

    void makeCoherentMapping( const unsigned int& nBins,
                              const std::function<std::vector<double>( const Event& )>& functors = nullptr,
                              const size_t& maxDepth = 20, const size_t& minPop = 50,
                              const std::function<uint64_t( const Event& )>& flagFunctor = nullptr );
    void makeCoherentMapping( const unsigned int& nBins, const BinDT& binDT );

    void groupByStrongPhase( const std::vector<CoherenceEvent>& coherenceEvents = {} );

    void readBinsFromFile( const std::string& binName ); ///
    unsigned int getBinNumber( const Event& event ) const;
    double operator()();
    void testCoherenceFactor() const;

    std::vector<HadronicParameters> getCoherenceFactors() const;
    std::vector<HadronicParameters> getCoherenceFactorsFast() const;
    std::vector<HadronicParameters> coherenceFactorsPerVoxel() const { return m_paramsPerVoxel; }

    void calculateNorms();
    void calculateCoherenceFactorsPerVoxel();
    void calculateCoherenceFactorsPerVoxel( const std::vector<CoherenceEvent>& dtEvents );

    void setGlobals( const double& globalPhase, const double& globalR )
    {
      m_globalR     = globalR;
      m_globalPhase = globalPhase;
    }
    void setBinOfNode( const unsigned int& nodeID, const unsigned int& bin ) { m_nodeID2Bin[nodeID] = bin; }
  };
} // namespace AmpGen

#endif
