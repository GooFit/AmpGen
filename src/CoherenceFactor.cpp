#include "AmpGen/CoherenceFactor.h"

#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/DynamicContainer.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/PhaseSpace.h"

using namespace AmpGen;

void CoherenceFactor::makeCoherentMapping( const unsigned int& nBins,
                                           const std::function<std::vector<double>( const Event& )>& functors,
                                           const size_t& maxDepth, const size_t& minPop,
                                           const std::function<uint64_t( const Event& )>& flagFunctor )
{
  if ( m_data != nullptr ) {
    m_voxels = BinDT( *m_data, Functor( functors ), MinEvents( minPop ) );
    return;
  }
  size_t population = pow( 2, maxDepth ) * 1.2 * minPop; /// 20 % safety factor ///
  Generator<> phsp( m_type );
  phsp.setRandom( new TRandom3() );
  size_t counter   = 0;
  int otherCounter = 0;
  size_t blockSize = 1 * pow( 10, 6 );
  auto makeEvents  = [&]( EventList& eventList ) {
    phsp.fillEventListPhaseSpace( eventList, blockSize, m_pdf1->size() + m_pdf2->size() );
    m_pdf1->reset( true );
    m_pdf2->reset( true );
    m_pdf1->setEvents( eventList );
    m_pdf2->setEvents( eventList );
    m_pdf1->prepare();
    m_pdf2->prepare();
    counter += blockSize;
    INFO( "Generated: " << counter << " DT events (out of " << population
                        << " requested); [evt counter=" << otherCounter << "]" );
  };

  DynamicContainer<Event, EventList> dlist( population, makeEvents );

  std::vector<CoherenceEvent> DTevents( population );
  std::vector<double*> addresses;
  std::vector<double*> addresses_veto;
  for ( auto evt : dlist ) {
    auto& te    = DTevents[otherCounter++];
    auto coords = functors( evt );
    for ( unsigned int i = 0; i < 5; ++i ) te.values[i] = coords[i];
    te.weight                                           = evt.weight() / evt.genPdf();
    te.amp1                                             = m_pdf1->getVal( evt );
    te.amp2                                             = m_pdf2->getVal( evt );
    if ( flagFunctor != nullptr && flagFunctor( evt ) ) {
      te.flag = true;
      addresses_veto.push_back( te.values.data() );
    } else {
      addresses.push_back( te.values.data() );
    }
  }
  m_voxels = BinDT( addresses, Functor( functors ), MaxDepth( maxDepth ), MinEvents( minPop ) );

  INFO( "Generated " << dlist.size()
                     << " events; entries per node =  " << double( dlist.size() ) / double( m_voxels.nodes().size() ) );
  m_nBins = nBins;
  groupByStrongPhase( DTevents );
}

void CoherenceFactor::makeCoherentMapping( const unsigned int& nBins, const BinDT& binDT )
{
  m_voxels = binDT;
  m_nBins  = nBins;
  groupByStrongPhase();
}

void CoherenceFactor::calculateCoherenceFactorsPerVoxel()
{
  m_paramsPerVoxel.clear();
  m_paramsPerVoxel.resize( m_voxels.size() + 1 );
  for ( unsigned int i = 0; i < m_paramsPerVoxel.size(); ++i ) m_paramsPerVoxel[i].id = i;
  auto pdf1CacheMap = m_pdf1->cacheAddresses( *m_data );
  auto pdf2CacheMap = m_pdf2->cacheAddresses( *m_data );

  for ( auto& evt : *m_data ) {
    int binNumber           = m_voxels.getBinNumber( evt );
    double w                = evt.weight() / evt.genPdf();
    std::complex<double> p1 = m_pdf1->getVal( evt, pdf1CacheMap );
    std::complex<double> p2 = m_pdf2->getVal( evt, pdf2CacheMap );
    m_paramsPerVoxel[binNumber].add( p1, p2, w );
  }
}

void CoherenceFactor::calculateCoherenceFactorsPerVoxel( const std::vector<CoherenceEvent>& dtEvents )
{
  m_paramsPerVoxel.clear();
  m_paramsPerVoxel.resize( m_voxels.size() + 1 );
  INFO( "h-Voxels = " << m_voxels.size() << " associating with ids" );
  for ( unsigned int i = 0; i < m_paramsPerVoxel.size(); ++i ) m_paramsPerVoxel[i].id = i;
  for ( auto& evt : dtEvents ) {
    unsigned int binNumber = evt.flag ? m_voxels.size() : m_voxels.getBinNumber( evt.values.data() );
    DEBUG( "BinNumber: " << binNumber << " pos = {" << evt.values[0] << "," << evt.values[1] << "," << evt.values[2]
                         << "," << evt.values[3] << "," << evt.values[4] << "}" );
    m_paramsPerVoxel[binNumber].add( evt.amp1, evt.amp2, evt.weight );
  }
}

void CoherenceFactor::groupByStrongPhase( const std::vector<CoherenceEvent>& coherenceEvents )
{

  if ( m_data != nullptr )
    calculateCoherenceFactorsPerVoxel();
  else
    calculateCoherenceFactorsPerVoxel( coherenceEvents );
  INFO( "Got " << m_voxels.size() << " voxels" );
  HadronicParameters global_hp =
      std::accumulate( m_paramsPerVoxel.begin(), m_paramsPerVoxel.end() - 1, HadronicParameters() );
  for ( auto& h : m_paramsPerVoxel ) global_hp += h;

  INFO( "Global parameters = " << global_hp.to_string() );
  double phiTransform = m_globalPhase - global_hp.d();
  double rTransform   = m_globalR * m_globalR * global_hp.n2 / global_hp.n1;

  INFO( "Global rotation = " << phiTransform * 180 / M_PI );
  for ( auto& h : m_paramsPerVoxel ) {
    h.scale_p1( 1. / rTransform );
    h.rotate( -phiTransform );
  }
  INFO( "Transform = " << sqrt( rTransform ) << " ; phi = " << phiTransform );

  HadronicParameters global_hp_new =
      std::accumulate( m_paramsPerVoxel.begin(), m_paramsPerVoxel.end() - 1, HadronicParameters() );
  INFO( global_hp_new.to_string() );
  std::sort( m_paramsPerVoxel.begin(), m_paramsPerVoxel.end() - 1,
             []( auto& b1, auto& b2 ) { return b1.d() > b2.d(); } );

  double norm             = global_hp_new.n1 / double( m_nBins );
  double currentNorm      = 0;
  unsigned int currentBin = 0;
  std::vector<std::pair<double, double>> binLimits;
  std::vector<unsigned int> voxels_in_this_bin;
  binLimits.emplace_back( -M_PI, 0 );
  voxels_in_this_bin.push_back( 0 );

  auto currentBinLimits = binLimits.rbegin();
  HadronicParameters hp_test;
  HadronicParameters glob_vox_test;
  for ( auto ihp = m_paramsPerVoxel.begin(); ihp != m_paramsPerVoxel.end() - 1; ++ihp ) {
    auto& hp = *ihp;
    currentNorm += hp.n1;
    hp.binID            = currentBin;
    m_nodeID2Bin[hp.id] = currentBin;
    ( *voxels_in_this_bin.rbegin() )++;
    hp_test += hp;
    double maxPhase = hp.d();
    if ( currentNorm > norm ) {
      currentNorm              = 0;
      currentBinLimits->second = maxPhase;
      currentBin++;
      binLimits.emplace_back( maxPhase, 0 );
      currentBinLimits = binLimits.rbegin();
      voxels_in_this_bin.push_back( 0 );
      INFO( "[ " << currentBin << " ]    " << hp_test.to_string() );
      glob_vox_test += hp_test;
      hp_test.reset();
    }
  }
  glob_vox_test += hp_test;
  m_paramsPerVoxel.rbegin()->binID = m_nBins; // allocate veto bin ///
  INFO( "[ " << currentBin + 1 << " ]    " << hp_test.to_string() );
  INFO( "[ VETO     ] " << m_paramsPerVoxel.rbegin()->to_string() );
  INFO( "[ NO VETO  ] " << glob_vox_test.to_string() );

  glob_vox_test += *m_paramsPerVoxel.rbegin();
  INFO( "[SUM VOXELS] " << global_hp_new.to_string() );
  INFO( "[SUM BINS  ] " << glob_vox_test.to_string() );

  for ( auto& limits : binLimits ) {
    INFO( "Bin Limits [" << limits.first << ", " << limits.second << "]" );
  }
}

void CoherenceFactor::testCoherenceFactor() const
{

  BinnedIntegrator<8> bid( m_data );
  bid.setView( [this]( auto& evt ) { return this->getBinNumber( evt ); } );
  std::vector<CoherenceCalculator> calcs( 8, CoherenceCalculator( m_pdf1, m_pdf2 ) );

  size_t PDF1_size = m_pdf1->size();
  size_t PDF2_size = m_pdf2->size();
  for ( unsigned int i = 0; i < PDF1_size; ++i ) {
    for ( unsigned int j = 0; j < PDF2_size; ++j ) {
      bid.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf2 )[j].pdf, [i, j, &calcs]( const auto& val ) {
        for ( unsigned int bin = 0; bin < 8; ++bin ) calcs[bin].R.set( i, j, val[bin] );
      } );
    }
  }
  for ( unsigned int i = 0; i < PDF1_size; ++i ) {
    for ( unsigned int j = i; j < PDF1_size; ++j ) {
      bid.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf1 )[j].pdf, [i, j, &calcs]( const auto& val ) {
        for ( unsigned int bin = 0; bin < 8; ++bin ) {
          calcs[bin].N1.set( i, j, val[bin] );
          if ( i != j ) calcs[bin].N1.set( j, i, std::conj( val[bin] ) );
        }
      } );
    }
  }
  for ( unsigned int i = 0; i < PDF2_size; ++i ) {
    for ( unsigned int j = i; j < PDF2_size; ++j ) {
      bid.addIntegral( ( *m_pdf2 )[i].pdf, ( *m_pdf2 )[j].pdf, [i, j, &calcs]( const auto& val ) {
        for ( unsigned int bin = 0; bin < 8; ++bin ) {
          calcs[bin].N2.set( i, j, val[bin] );
          if ( i != j ) calcs[bin].N2.set( j, i, std::conj( val[bin] ) );
        }
      } );
    }
  }
  bid.flush();

  for ( unsigned int i = 0; i < 8; ++i ) {
    std::complex<double> VALUE = calcs[i].getVal();
    INFO( "Bin [" << i + 1 << "] R=" << std::abs( VALUE ) << " d= " << 180 * std::arg( VALUE ) / M_PI
                  << " N1= " << calcs[i].getN1() << " N2= " << calcs[i].getN2() );
  }
}

std::vector<HadronicParameters> CoherenceFactor::getCoherenceFactors() const
{
  std::vector<HadronicParameters> R( m_nBins + 1 );
  HadronicParameters RG;
  auto pdf1CacheMap = m_pdf1->cacheAddresses( *m_data );
  auto pdf2CacheMap = m_pdf2->cacheAddresses( *m_data );
  for ( auto& evt : *m_data ) {
    unsigned int binNumber  = getBinNumber( evt );
    double w                = evt.weight() / evt.genPdf();
    std::complex<double> p1 = m_pdf1->getVal( evt, pdf1CacheMap );
    std::complex<double> p2 = m_pdf2->getVal( evt, pdf2CacheMap );
    R[binNumber].add( p1, p2, w );
    RG.add( p1, p2, w );
  }
  double phiTransform = m_globalPhase - std::arg( RG.coherence );
  double rTransform   = m_globalR * m_globalR * RG.n2 / RG.n1;
  for ( auto& h : R ) {
    h.scale_p1( 1. / rTransform );
    h.rotate( -phiTransform );
  }
  INFO( RG.to_string() );
  return R;
}

std::vector<HadronicParameters> CoherenceFactor::getCoherenceFactorsFast() const
{

  std::vector<HadronicParameters> coherence_factors( m_nBins + 1 );
  for ( auto hp = m_paramsPerVoxel.begin(); hp != m_paramsPerVoxel.end() - 1; ++hp ) {
    coherence_factors[hp->binID] += *hp;
  }
  coherence_factors[m_nBins] = *m_paramsPerVoxel.rbegin();
  return coherence_factors;
}

std::vector<size_t> CoherenceFactor::getNumberOfEventsInEachBin( const EventList& events ) const
{

  std::vector<size_t> counter( m_nBins + 1, 0 );
  for ( auto& event : events ) {
    counter[getBinNumber( event )]++;
  }
  return counter;
}

CoherenceFactor::CoherenceFactor( const std::string& filename ) { readBinsFromFile( filename ); }

void CoherenceFactor::readBinsFromFile( const std::string& binName )
{
  INFO( "Reading file = " << binName );
  m_voxels = BinDT( File( binName ) );
  m_nBins  = 0;
  for ( auto& node : m_voxels.nodes() ) {
    m_nodeID2Bin[node->voxNumber()]                                         = node->binNumber();
    if ( node->binNumber() != 999 && node->binNumber() >= m_nBins ) m_nBins = node->binNumber() + 1;
  }
  INFO( "Got " << m_nBins << " coherent bins and " << m_voxels.size() << " voxels" );
}

double CoherenceFactor::phase() { return std::arg( m_global.getCoherence() ); }
double CoherenceFactor::getR() { return std::abs( m_global.getCoherence() ); }

HadronicParameters CoherenceFactor::calculateGlobalCoherence( EventList* _data )
{
  m_global.clear();
  bool freshData = false;
  if ( _data == nullptr ) {
    _data = new EventList( m_type );
    PhaseSpace phsp( m_type );
    phsp.setRandom( new TRandom3() );
    for ( unsigned int i = 0; i < 1e6; ++i ) _data->push_back( phsp.makeEvent( m_pdf1->size() + m_pdf2->size() ) );
    freshData            = true;
  }

  m_pdf1->setMC( *_data );
  m_pdf2->setMC( *_data );
  m_pdf1->prepare();
  m_pdf2->prepare();
  auto pdf1CacheMap = m_pdf1->cacheAddresses( *_data );
  auto pdf2CacheMap = m_pdf2->cacheAddresses( *_data );
  for ( unsigned int i = 0; i < _data->size(); ++i ) {
    const Event& evt          = ( *_data )[i];
    double w                  = evt.weight() / evt.genPdf();
    std::complex<double> pdf1 = m_pdf1->getVal( evt, pdf1CacheMap );
    std::complex<double> pdf2 = m_pdf2->getVal( evt, pdf2CacheMap );
    m_global.add( pdf1, pdf2, w );
  }

  if ( freshData ) delete _data;
  return m_global;
}

void CoherenceFactor::writeToFile( const std::string& filename )
{
  auto endNodes = m_voxels.nodes();
  for ( auto& node : endNodes ) {
    node->setBinNumber( m_nodeID2Bin[node->voxNumber()] );
  }
  m_voxels.serialize( filename );
}

HadronicParameters CoherenceFactor::getVal() const { return m_global; }
CoherenceFactor::CoherenceFactor() = default;

CoherenceFactor::CoherenceFactor( FastCoherentSum* pdf1, FastCoherentSum* pdf2, EventList* data )
    : m_pdf1( pdf1 ), m_pdf2( pdf2 ), m_calc( pdf1, pdf2 ), m_data( data )
{
  if ( m_data != nullptr ) {
    calculateGlobalCoherence( m_data );
    calculateNorms();
    m_type = m_data->eventType();
  }
}

CoherenceFactor::CoherenceFactor( FastCoherentSum* pdf1, FastCoherentSum* pdf2, const EventType& type )
    : m_pdf1( pdf1 ), m_pdf2( pdf2 ), m_calc( pdf1, pdf2 ), m_data( nullptr ), m_type( type )
{
}

void CoherenceFactor::calculateNorms()
{
  Integrator<10> id( m_data );
  for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
    for ( unsigned int j = 0; j < m_pdf2->size(); ++j ) {
      id.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf2 )[j].pdf,
                      [i, j, this]( const std::complex<double>& val ) { this->m_calc.R.set( i, j, val ); } );
    }
  }

  for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
    for ( unsigned int j = i; j < m_pdf1->size(); ++j ) {
      id.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf1 )[j].pdf, [i, j, this]( const std::complex<double>& val ) {
        this->m_calc.N1.set( i, j, val );
        if ( i != j ) this->m_calc.N1.set( j, i, std::conj( val ) );
      } );
    }
  }

  for ( unsigned int i = 0; i < m_pdf2->size(); ++i ) {
    for ( unsigned int j = i; j < m_pdf2->size(); ++j ) {
      id.addIntegral( ( *m_pdf2 )[i].pdf, ( *m_pdf2 )[j].pdf, [i, j, this]( const std::complex<double>& val ) {
        this->m_calc.N2.set( i, j, val );
        if ( i != j ) this->m_calc.N2.set( j, i, std::conj( val ) );
      } );
    }
  }

  id.flush();
}

unsigned int CoherenceFactor::getBinNumber( const Event& event ) const
{
  int voxelID = m_voxels.getBinNumber( event );
  return voxelID == -1 ? m_nBins : m_nodeID2Bin.find( voxelID )->second;
}
double CoherenceFactor::operator()() { return std::abs( m_calc.getVal() ); }
