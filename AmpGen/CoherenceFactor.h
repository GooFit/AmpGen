#ifndef AMPGEN_COHERENCEFACTOR_H
#define AMPGEN_COHERENCEFACTOR_H
#include <bits/stdint-uintn.h>
#include <stddef.h>
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
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class Event;

  struct CoherenceEvent {
    std::array<real_t, 5> values;
    complex_t amp1;
    complex_t amp2;
    real_t weight;
    uint64_t flag;
  };

  struct CoherenceCalculator {
    CoherentSum* m_pdf1;
    CoherentSum* m_pdf2;
    Bilinears R;
    Bilinears N1;
    Bilinears N2;
    real_t rTransform;
    real_t phiTransform;
    CoherenceCalculator() = default;
    CoherenceCalculator( CoherentSum* pdf1, CoherentSum* pdf2 )
        : m_pdf1( pdf1 )
        , m_pdf2( pdf2 )
        , R( m_pdf1->size(), m_pdf2->size() )
        , N1( m_pdf1->size(), m_pdf1->size() )
        , N2( m_pdf2->size(), m_pdf2->size() )
        , rTransform( 1 )
        , phiTransform( 0 )
    {
    }

    void setTransform( const real_t& rT, const real_t& pT )
    {
      rTransform   = rT;
      phiTransform = pT;
    }
    real_t getN2() const
    {
      real_t accN2 = 0;
      for ( unsigned int i = 0; i < m_pdf2->size(); ++i ) {
        accN2 += std::norm( ( *m_pdf2 )[i].coefficient ) * std::abs( N2.get( i, i ) );
        for ( unsigned int j = i + 1; j < m_pdf2->size(); ++j ) {
          accN2 +=
              2 * std::real( ( *m_pdf2 )[i].coefficient * std::conj( ( *m_pdf2 )[j].coefficient ) * N2.get( i, j ) );
        }
      }
      return accN2;
    }

    real_t getN1() const
    {
      real_t accN1 = 0;
      for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
        accN1 += std::norm( ( *m_pdf1 )[i].coefficient ) * std::abs( N1.get( i, i ) );
        for ( unsigned int j = i + 1; j < m_pdf1->size(); ++j ) {
          accN1 +=
              2 * std::real( ( *m_pdf1 )[i].coefficient * std::conj( ( *m_pdf1 )[j].coefficient ) * N1.get( i, j ) );
        }
      }
      return accN1 * rTransform * rTransform;
    }

    complex_t getVal()
    {
      complex_t accR( 0, 0 );
      m_pdf1->transferParameters();
      m_pdf2->transferParameters();
      for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
        for ( unsigned int j = 0; j < m_pdf2->size(); ++j ) {
          accR += ( *m_pdf1 )[i].coefficient * std::conj( ( *m_pdf2 )[j].coefficient ) * R.get( i, j );
        }
      }
      return complex_t( cos( phiTransform ), sin( phiTransform ) ) * rTransform * accR /
             sqrt( getN1() * getN2() );
    }
  };

  struct HadronicParameters {

    HadronicParameters( const real_t& _R, const real_t& _d, const real_t& k1 );
    HadronicParameters( const real_t& _R, const real_t& _d, const real_t& k1, const real_t& k2 );
    HadronicParameters() : n1( 0 ), n2( 0 ), wt(0), coherence( 0, 0 ), nFills( 0 ), id( 0 ), binID( 0 ) {}

    real_t n1;
    real_t n2;
    real_t wt;
    complex_t coherence;
    unsigned int nFills;
    unsigned int id;
    unsigned int binID;
    real_t I1() const { return n1 / wt  ; }
    real_t I2() const { return n2 / wt  ; }
    real_t d() const { return std::arg( coherence ); }
    real_t R() const { return std::abs( getCoherence() ); }
    real_t r() const { return sqrt( n1 / n2 ); } 
    complex_t getCoherence() const { return coherence / ( wt * sqrt( I1() * I2() ) ) ; }
    void scale_p1( const real_t& sf );
    void scale_p2( const real_t& sf );
    void rotate( const real_t& angle );
    void clear();

    void add( const complex_t& f1, const complex_t& f2, const real_t& weight = 1 );

    HadronicParameters operator+=( const HadronicParameters& other );
    HadronicParameters operator-=( const HadronicParameters& other );

    std::string to_string() const;
  };

  HadronicParameters operator+( const HadronicParameters& a, const HadronicParameters& b );
  class CoherenceFactor
  {
  private:
    HadronicParameters m_global;
    real_t m_globalPhase;
    real_t m_globalR;
    CoherentSum* m_pdf1;
    CoherentSum* m_pdf2;
    BinDT m_voxels;
    std::vector<HadronicParameters> m_paramsPerVoxel;
    std::map<unsigned int, unsigned int> m_nodeID2Bin;
    unsigned int m_nBins;

    CoherenceCalculator m_calc;
    EventList* m_data;
    EventType m_type;

    void MakeEventDeltaPlots(const std::vector<CoherenceEvent>& events );
  public:
    BinDT& getDT() { return m_voxels; }
    real_t phase();
    real_t getR();
    std::vector<size_t> getNumberOfEventsInEachBin( const EventList& events ) const;
    HadronicParameters calculateGlobalCoherence( EventList* data = nullptr );
    void writeToFile( const std::string& filename );

    HadronicParameters getVal() const;
    CoherenceFactor( const std::string& filename );
    CoherenceFactor();
    CoherenceFactor( CoherentSum* pdf1, CoherentSum* pdf2, EventList* data = nullptr );

    CoherenceFactor( CoherentSum* pdf1, CoherentSum* pdf2, const EventType& type );

    void makeCoherentMapping( const unsigned int& nBins,
                              const std::function<std::vector<real_t>( const Event& )>& functors = nullptr,
                              const size_t& maxDepth = 20, const size_t& minPop = 50,
                              const std::function<uint64_t( const Event& )>& flagFunctor = nullptr,
                              const bool& refreshEvents = true );
    void makeCoherentMapping( const unsigned int& nBins, const BinDT& binDT );

    void groupByStrongPhase( const std::vector<CoherenceEvent>& coherenceEvents = {}, const bool& recalculateVoxels =true );

    void readBinsFromFile( const std::string& binName ); ///
    unsigned int getBinNumber( const Event& event ) const;
    real_t operator()();
    void testCoherenceFactor() const;

    std::vector<HadronicParameters> getCoherenceFactors() const;
    std::vector<HadronicParameters> getCoherenceFactorsFast() const;
    std::vector<HadronicParameters> coherenceFactorsPerVoxel() const { return m_paramsPerVoxel; }

    void calculateNorms();
    void calculateCoherenceFactorsPerVoxel();
    void calculateCoherenceFactorsPerVoxel( const std::vector<CoherenceEvent>& dtEvents );

    void setGlobals( const real_t& globalPhase, const real_t& globalR )
    {
      m_globalR     = globalR;
      m_globalPhase = globalPhase;
    }
    void setBinOfNode( const unsigned int& nodeID, const unsigned int& bin ) { m_nodeID2Bin[nodeID] = bin; }
  };
} // namespace AmpGen

#endif
