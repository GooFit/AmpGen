#include "AmpGen/CoherenceFactor.h"

#include <TRandom3.h>
#include <memory.h>
#include <ext/alloc_traits.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>
#include <cstdint>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/DynamicContainer.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Event.h"
#include "AmpGen/Utilities.h"
#include "TH1.h"

using namespace AmpGen;

void HadronicParameters::scale_p1( const real_t& sf )
{
  n1 /= sf;
  coherence /= sqrt( sf );
}
void HadronicParameters::scale_p2( const real_t& sf )
{
  n2 /= sf;
  coherence /= sqrt( sf );
}

void HadronicParameters::rotate( const real_t& angle ) 
{ 
  coherence *= complex_t( cos( angle ), sin( angle ) ); 
}

void HadronicParameters::clear() 
{ 
  n1        = 0;
  n2        = 0;
  nFills    = 0;
  id        = 0;
  wt        = 0;
  coherence = complex_t( 0, 0 ); 
}

void HadronicParameters::add( const complex_t& f1, 
    const complex_t& f2, 
    const real_t& weight )
{
  n1        += weight * std::norm( f1 );
  n2        += weight * std::norm( f2 );
  wt        += weight; 
  coherence += weight * std::conj(f1) * f2 ;
  nFills++;
}

HadronicParameters::HadronicParameters( const real_t& _R, const real_t& _d, const real_t& k1 ) : 
  n1(k1*k1), 
  n2(1), 
  wt(1), 
  coherence(k1*_R*cos( M_PI*_d/180.), k1*_R*sin(M_PI*_d/180.)) {}

HadronicParameters::HadronicParameters( const real_t& _R, const real_t& _d, const real_t& k1, const real_t& k2 )
  : n1( k1 ), 
  n2( k2 ), 
  wt(1), 
  coherence( sqrt( k1 * k2 ) * _R * cos( _d ), sqrt( k1 * k2 ) * _R * sin( _d ) ) {}

HadronicParameters HadronicParameters::operator+=( const HadronicParameters& other )
{
  n1        += other.n1;
  n2        += other.n2;
  coherence += other.coherence;
  nFills    += other.nFills;
  wt        += other.wt;
  return *this;
}

HadronicParameters HadronicParameters::operator-=( const HadronicParameters& other )
{
  n1        -= other.n1;
  n2        -= other.n2;
  coherence -= other.coherence;
  nFills    -= other.nFills;
  wt        -= other.wt;
  return *this;
}
HadronicParameters AmpGen::operator+( const HadronicParameters& a, const HadronicParameters& b )
{
  HadronicParameters hp;
  hp.wt        = a.wt + b.wt;
  hp.n1        = a.n1 + b.n1 ;
  hp.n2        = a.n2 + b.n2 ;
  hp.coherence = a.coherence + b.coherence ;
  hp.nFills    = a.nFills + b.nFills;
  return hp;
}

std::string HadronicParameters::to_string() const
{
  auto c = getCoherence();
  std::string total = "R = " + round( R() , 3 ) + " δ = " + round( 180 * d() / M_PI, 2 ) +
    " r = " + round( r(), 4 ) + " K = " + round( I1() , 6 ) + " K' = " + round( I2(), 6 ) + " N = " +
    std::to_string( wt ) + " ci = " + std::to_string( std::real( c ) ) + " si = " +
    std::to_string( std::imag( c ) ) ; 
  return total;
}

void CoherenceFactor::makeCoherentMapping( const unsigned int& nBins,
    const std::function<std::vector<real_t>( const Event& )>& functors,
    const size_t& maxDepth, 
    const size_t& minPop,
    const std::function<uint64_t( const Event& )>& flagFunctor,
    const bool& refreshEvents )
{
  m_nBins = nBins;

  if ( m_data != nullptr ) {
    m_voxels = BinDT( *m_data, Functor( functors ), MinEvents( minPop ) );
    return;
  }
  if( !refreshEvents ){
    groupByStrongPhase({},false); 
    return;
  }

  size_t population = pow( 2, maxDepth ) * 1.2 * minPop; /// 20 % safety factor ///
  Generator<> phsp( m_type );
  phsp.setRandom( new TRandom3() );
  size_t counter   = 0;
  int otherCounter = 0;
  size_t blockSize = 1 * pow( 10, 5 );
  auto makeEvents  = [&]( EventList& eventList ) {
    bool useFlatEvents = NamedParameter<bool> ("CoherenceFactor::FlatEvents",true);
    if( useFlatEvents ) phsp.fillEventListPhaseSpace( eventList, blockSize, 0);
    else phsp.fillEventList( *m_pdf2, eventList, blockSize );

    eventList.resetCache();

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
  std::vector<real_t*> addresses;
  std::vector<real_t*> addresses_veto;
  for ( auto evt : dlist ) {
    auto& te    = DTevents[otherCounter++];
    te.weight   = 1./evt.genPdf();
    te.amp1     = m_pdf1->getVal( evt );
    te.amp2     = m_pdf2->getVal( evt );
    auto coords = functors( evt );
    for ( unsigned int i = 0; i < 5; ++i ) 
      te.values[i] = coords[i];
    if ( flagFunctor != nullptr && flagFunctor( evt ) ) {
      te.flag = true;
      addresses_veto.push_back( te.values.data() );
    } else {
      addresses.push_back( te.values.data() );
    }
  }
  m_voxels = BinDT( addresses, Functor( functors ), MaxDepth( maxDepth ), MinEvents( minPop ) );

  INFO( "Generated " << dlist.size()
      << " events; entries per node =  " << real_t( dlist.size() ) / real_t( m_voxels.nodes().size() ) );
  groupByStrongPhase( DTevents , true );
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

  TH1D* deltaPhaseWeightPdf1 = new TH1D( "deltaPhi1","",100,-180,180);
  TH1D* deltaPhaseWeightPdf2 = new TH1D( "deltaPhi2","",100,-180,180);
  TH1D* deltaPhaseNoWeight   = new TH1D( "deltaPhase0","",100,-180,180);
  for ( auto& evt : *m_data ) {
    int binNumber           = m_voxels.getBinNumber( evt );
    real_t w                = evt.weight() / evt.genPdf();
    complex_t p1 = m_pdf1->getVal( evt, pdf1CacheMap );
    complex_t p2 = m_pdf2->getVal( evt, pdf2CacheMap );
    real_t arg = std::arg( p1* std::conj(p2 ) );
    deltaPhaseWeightPdf1->Fill( arg * 180 / M_PI, std::abs(p1) );
    deltaPhaseWeightPdf2->Fill( arg * 180 / M_PI, std::abs(p2) );
    deltaPhaseNoWeight->Fill( arg * 180 / M_PI, 1 );
    m_paramsPerVoxel[binNumber].add( p1, p2, w );
  }
  deltaPhaseWeightPdf1->Write();
  deltaPhaseWeightPdf2->Write();
  deltaPhaseNoWeight->Write();
}

void CoherenceFactor::calculateCoherenceFactorsPerVoxel( const std::vector<CoherenceEvent>& dtEvents )
{
  m_paramsPerVoxel.clear();
  m_paramsPerVoxel.resize( m_voxels.size() + 1 );
  INFO( "h-Voxels = " << m_voxels.size() << " associating with ids" );
  for ( unsigned int i = 0; i < m_paramsPerVoxel.size(); ++i ) m_paramsPerVoxel[i].id = i;

  TH1D* deltaPhaseNoWeight  = new TH1D( "deltaPhi0","",100,-180,180);
  TH1D* deltaPhaseWeightVN = new TH1D( "deltaVPhiN","",100,-180,180);
  TH1D* weight = new TH1D("weight","",100,0,100);

  HadronicParameters hG; 
  for ( auto& evt : dtEvents ) {
    unsigned int binNumber = evt.flag ? m_voxels.size() : m_voxels.getBinNumber( evt.values.data() );
    DEBUG( "BinNumber: " << binNumber << " pos = {" << evt.values[0] << "," << evt.values[1] << "," << evt.values[2]
        << "," << evt.values[3] << "," << evt.values[4] << "}" );
    real_t arg = std::arg( std::conj(evt.amp1) * evt.amp2 );
    deltaPhaseNoWeight->Fill(  arg * 180 / M_PI, 1 );
    weight->Fill( evt.weight );
    m_paramsPerVoxel[binNumber].add( evt.amp1, evt.amp2, evt.weight );
    hG.add( evt.amp1, evt.amp2, evt.weight );
  }
  INFO( hG.to_string() );
  for( auto& v : m_paramsPerVoxel ) deltaPhaseWeightVN->Fill( v.d()*180/M_PI, 1 );
  deltaPhaseNoWeight->Write();
  deltaPhaseWeightVN->Write();
  weight->Write();
}

void CoherenceFactor::groupByStrongPhase( const std::vector<CoherenceEvent>& coherenceEvent, const bool& recalculateVoxels )
{
  if( recalculateVoxels ){
    if ( m_data != nullptr )
      calculateCoherenceFactorsPerVoxel();
    else
      calculateCoherenceFactorsPerVoxel( coherenceEvent );
  }
  INFO( "Got " << m_voxels.size() << " voxels" );
  HadronicParameters global_hp;
  for ( auto& h : m_paramsPerVoxel ) global_hp += h;

  INFO( "Global parameters = " << global_hp.to_string() );
  real_t phiTransform = m_globalPhase - global_hp.d();
  real_t rTransform   = m_globalR * m_globalR * global_hp.n2 / global_hp.n1;

  INFO( "Global rotation = " << phiTransform * 180 / M_PI );
  for ( auto& h : m_paramsPerVoxel ) {
    h.scale_p1( 1. / rTransform );
    h.rotate( phiTransform );
  }
  INFO( "Transform = " << sqrt( rTransform ) << " ; phi = " << phiTransform );

  HadronicParameters global_hp_new = std::accumulate( m_paramsPerVoxel.begin(), m_paramsPerVoxel.end() - 1, HadronicParameters() );
  INFO( global_hp_new.to_string() );
  std::sort( m_paramsPerVoxel.begin(), m_paramsPerVoxel.end() - 1, []( auto& b1, auto& b2 ) { return b1.d() < b2.d(); } );

  real_t norm_target      = global_hp_new.wt / real_t( m_nBins );
  unsigned int currentBin = 0;
  std::vector<std::pair<real_t, real_t>> binLimits;
  binLimits.emplace_back( -M_PI, 0 );
  INFO( "Normalisation-per-bin = " << norm_target << " i1 = " << global_hp_new.I1()  );
  auto currentBinLimits = binLimits.rbegin();
  HadronicParameters hp_accumulator, glob_vox_test;

  for( auto ihp = m_paramsPerVoxel.begin(); ihp != m_paramsPerVoxel.end() - 1; ++ihp ) {
    ihp->binID            = currentBin;
    m_nodeID2Bin[ihp->id] = currentBin;
    hp_accumulator += *ihp;
    if ( hp_accumulator.wt > norm_target ) {
      currentBinLimits->second = ihp->d();
      currentBin++;
      binLimits.emplace_back( ihp->d(), 0  );
      currentBinLimits = binLimits.rbegin();
      glob_vox_test += hp_accumulator;
      hp_accumulator.clear();
    }
  }
  INFO("Generated: " << currentBin << " bins");
  glob_vox_test += hp_accumulator;
  m_paramsPerVoxel.rbegin()->binID = m_nBins; // allocate veto bin ///
  //  for( auto& bin : m_para
  // INFO( "[ " << currentBin + 1 << "       ]    " << hp_test.to_string() );
  INFO( "[ VETO     ] " << m_paramsPerVoxel.rbegin()->to_string() );
  INFO( "[ NO VETO  ] " << glob_vox_test.to_string() );

  glob_vox_test += *m_paramsPerVoxel.rbegin();
  INFO( "[SUM VOXELS] " << global_hp_new.to_string() );
  INFO( "[SUM BINS  ] " << glob_vox_test.to_string() );
  if( recalculateVoxels  && m_data == nullptr ) MakeEventDeltaPlots( coherenceEvent );
  //  for ( auto& limits : binLimits ) {
  //    INFO( "Bin Limits [" << limits.first << ", " << limits.second << "]" );
  //  }
}

void CoherenceFactor::MakeEventDeltaPlots( const std::vector<CoherenceEvent>& events ) {

  std::vector<TH1D*> plots; 
  for(size_t i = 0 ; i < m_nBins ; ++i ) 
    plots.push_back( new TH1D( ("Bin_"+std::to_string(i)+"_deltaPhi").c_str(),0,100,-180,180) );
  for( auto& event : events ){
    unsigned int voxNumber = event.flag ? m_voxels.size() : m_voxels.getBinNumber( event.values.data() );
    real_t arg = std::arg( std::conj(event.amp1) * event.amp2 );
    auto binNumber = m_nodeID2Bin[ voxNumber ];
    if( ! event.flag ) plots[binNumber]->Fill( arg*180/M_PI );
  }
  for( auto& plot : plots ) plot->Write();
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
    complex_t VALUE = calcs[i].getVal();
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
    real_t w                = evt.weight() / evt.genPdf();
    complex_t p1 = m_pdf1->getVal( evt, pdf1CacheMap );
    complex_t p2 = m_pdf2->getVal( evt, pdf2CacheMap );
    R[binNumber].add( p1, p2, w );
    RG.add( p1, p2, w );
  }
  real_t phiTransform = m_globalPhase - std::arg( RG.coherence );
  real_t rTransform   = m_globalR * m_globalR * RG.n2 / RG.n1;
  for ( auto& h : R ) {
    h.scale_p1( 1. / rTransform );
    h.rotate(  phiTransform );
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

real_t CoherenceFactor::phase() { return std::arg( m_global.getCoherence() ); }
real_t CoherenceFactor::getR() { return std::abs( m_global.getCoherence() ); }

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
    real_t w                  = evt.weight() / evt.genPdf();
    complex_t pdf1 = m_pdf1->getVal( evt, pdf1CacheMap );
    complex_t pdf2 = m_pdf2->getVal( evt, pdf2CacheMap );
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

CoherenceFactor::CoherenceFactor( CoherentSum* pdf1, CoherentSum* pdf2, EventList* data )
  : m_pdf1( pdf1 ), m_pdf2( pdf2 ), m_calc( pdf1, pdf2 ), m_data( data )
{
  if ( m_data != nullptr ) {
    calculateGlobalCoherence( m_data );
    calculateNorms();
    m_type = m_data->eventType();
  }
}

CoherenceFactor::CoherenceFactor( CoherentSum* pdf1, CoherentSum* pdf2, const EventType& type )
  : m_pdf1( pdf1 ), m_pdf2( pdf2 ), m_calc( pdf1, pdf2 ), m_data( nullptr ), m_type( type )
{
}

void CoherenceFactor::calculateNorms()
{
  Integrator<10> id( m_data );
  for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
    for ( unsigned int j = 0; j < m_pdf2->size(); ++j ) {
      id.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf2 )[j].pdf,
          [i, j, this]( const complex_t& val ) { this->m_calc.R.set( i, j, val ); } );
    }
  }

  for ( unsigned int i = 0; i < m_pdf1->size(); ++i ) {
    for ( unsigned int j = i; j < m_pdf1->size(); ++j ) {
      id.addIntegral( ( *m_pdf1 )[i].pdf, ( *m_pdf1 )[j].pdf, [i, j, this]( const complex_t& val ) {
          this->m_calc.N1.set( i, j, val );
          if ( i != j ) this->m_calc.N1.set( j, i, std::conj( val ) );
          } );
    }
  }

  for ( unsigned int i = 0; i < m_pdf2->size(); ++i ) {
    for ( unsigned int j = i; j < m_pdf2->size(); ++j ) {
      id.addIntegral( ( *m_pdf2 )[i].pdf, ( *m_pdf2 )[j].pdf, [i, j, this]( const complex_t& val ) {
          this->m_calc.N2.set( i, j, val );
          if ( i != j ) this->m_calc.N2.set( j, i, std::conj( val ) );
          } );
    }
  }

  id.flush();
}

std::vector<unsigned int> CoherentSum::cacheAddresses( const EventList& evts ) const
{
  std::vector<unsigned int> addresses;
  for ( auto& mE : m_matrixElements ) {
    addresses.push_back( evts.getCacheIndex( mE.pdf ) );
  }
  return addresses;
}

complex_t CoherentSum::getVal( const Event& evt ) const
{
  complex_t value( 0., 0. );
  for ( auto& mE : m_matrixElements ) {
    value += mE.coefficient * evt.getCache( mE.addressData );
  }
  return value;
}

complex_t CoherentSum::getVal( const Event& evt, const std::vector<unsigned int>& cacheAddresses ) const
{
  complex_t value( 0., 0. );
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i )
    value += m_matrixElements[i].coefficient * evt.getCache( cacheAddresses[i] );
  return value;
}

unsigned int CoherenceFactor::getBinNumber( const Event& event ) const
{
  int voxelID = m_voxels.getBinNumber( event );
  return voxelID == -1 ? m_nBins : m_nodeID2Bin.find( voxelID )->second;
}
real_t CoherenceFactor::operator()() { return std::abs( m_calc.getVal() ); }
