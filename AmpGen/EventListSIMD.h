#ifndef AMPGEN_EVENTLIST2_H
#define AMPGEN_EVENTLIST2_H

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/EventList.h"
#include <chrono>
#include <functional>
#include <numeric>
#include <cstddef>
#include <algorithm>

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#if ENABLE_AVX2

  #include "AmpGen/simd/avx2_types.h"
  #include "AmpGen/simd/iterator.h"

namespace AmpGen
{
  using float_v   = AVX2::float_t;
  using complex_v = AVX2::complex_t; 

  class CompiledExpressionBase; 
  class EventListSIMD
  {
  private:
    std::vector<float_v>         m_data              = {};
    std::vector<float_v>         m_weights           = {};
    std::vector<float_v>         m_genPDF            = {};
    std::vector<complex_v>       m_cache             = {};
    EventType                    m_eventType         = {};
    std::map<uint64_t, unsigned> m_pdfIndex          = {};
    unsigned                     m_eventSize         = {0};
    unsigned                     m_nEvents           = {0};
    unsigned                     m_nBlocks           = {0};
  public:
    typedef Event value_type;
    EventListSIMD() = default;
    EventListSIMD( const EventType& type );
    template < class ... ARGS > EventListSIMD( const std::string& fname, const EventType& evtType, const ARGS&... args ) : EventListSIMD(evtType) 
    {
      loadFromFile( fname, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventListSIMD( const std::string& fname, const ARGS&... args ) : EventListSIMD() 
    {
      loadFromFile( fname, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventListSIMD( const std::vector<std::string>& fname, const EventType& evtType, const ARGS&... args ) : EventListSIMD(evtType) 
    {
      for( auto& f : fname ) loadFromFile( f, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventListSIMD( TTree* tree, const EventType& evtType, const ARGS&... args ) : EventListSIMD(evtType)
    {
      loadFromTree( tree, ArgumentPack(args...) );
    }
    const float_v* data() const { return m_data.data(); }
    const AVX2::complex_t* cache() const { return m_cache.data() ; } 
    EventListSIMD( const EventList& other );     
    void resetCache();
    const AVX2::complex_t cache( const unsigned& evtIndex, const unsigned& cachePos )
    {
      return m_cache[ (unsigned)(evtIndex/float_v::size) * cacheSize() + cachePos ]; 
    }
    const Event at(const unsigned& p) const { return EventListSIMD::operator[](p) ; }
    const float_v* block(const unsigned& p) { return m_data.data() + p * m_eventSize ; }
    float_v weight(const unsigned& p) const { return m_weights[p]; }
    float_v genPDF(const unsigned& p) const { return m_genPDF[p]; }
    const Event operator[]( const size_t&) const;
    std::array<Event, float_v::size> scatter(unsigned) const;
    void gather(const std::array<Event, float_v::size>&, unsigned);   
    auto begin() const { return make_scatter_iterator<float_v::size>(0,this); }
    auto   end() const { return make_scatter_iterator<float_v::size>(m_nEvents, (const EventListSIMD*)(nullptr) ); } 
    auto begin()       { return make_scatter_iterator<float_v::size, true>(0, this); }
    auto   end()       { return make_scatter_iterator<float_v::size, true>(m_nEvents, (EventListSIMD*)(nullptr) ); }
    EventType eventType()                         const { return m_eventType; }
    size_t aligned_size()                         const { return nBlocks() * float_v::size; } ///aligned number of events
    size_t cacheSize()                            const { return m_cache.size() / m_nBlocks; }  /// number of cached elements
    double integral()                             const;
    size_t eventSize()                            const { return m_eventSize; }
    size_t size()                                 const { return m_nEvents ; }
    size_t nBlocks()                              const { return m_nBlocks; }
    void reserve( const size_t& size ) { m_data.reserve( size * m_eventType.size() ); }
    void setEventType( const EventType& type ) { m_eventType = type; m_eventSize = m_eventType.size(); }
    void add( const EventListSIMD& evts );
    void loadFromTree( TTree* tree, const ArgumentPack& args ); 
    void loadFromFile( const std::string& fname, const ArgumentPack& args );
    void printCacheInfo( const unsigned int& nEvt = 0 );
    void clear();

    TTree* tree( const std::string& name, const std::vector<std::string>& extraBranches = {} ) const;
    
    size_t getCacheIndex( const CompiledExpressionBase& PDF, bool& status ) const;
    size_t getCacheIndex( const CompiledExpressionBase& PDF ) const;
    template <class T> unsigned registerExpression(const T& expression, const unsigned& size_of=0)
    {
      auto key = FNV1a_hash( expression.name() );
      auto pdfIndex = m_pdfIndex.find( key );
      if ( pdfIndex != m_pdfIndex.end() ) return pdfIndex->second;
      else {
        unsigned nEvents = aligned_size();
        unsigned expression_size = size_of == 0 ? expression.returnTypeSize() / sizeof(AmpGen::AVX2::complex_t) : size_of; 
        m_pdfIndex[key] = m_cache.size() / nBlocks(); 
        m_cache.resize(m_cache.size() + nBlocks() * expression_size);
        return m_pdfIndex[key];
      }
    }
    template <class FCN> void updateCache( const FCN& fcn, const size_t& index )
    {
      fcn.batch(m_cache.data() + index, aligned_size(), m_eventSize, cacheSize(), fcn.externBuffer().data(), m_data.data());
    }
    void reserveCache(const unsigned& index);
    void resizeCache( const unsigned& index);
    TH1D* makeProjection(const Projection& projection  , const ArgumentPack& args = ArgumentPack()) const; 
    TH2D* makeProjection(const Projection2D& projection, const ArgumentPack& args = ArgumentPack()) const;
    std::vector<TH1D*> makeProjections( const std::vector<Projection>& projections, const ArgumentPack& args );

    template <class... ARGS> std::vector<TH1D*> makeDefaultProjections( const ARGS&... args )
    {
      auto argPack = ArgumentPack( args... );
      size_t nBins = argPack.getArg<PlotOptions::Bins>(100);
      auto proj = eventType().defaultProjections(nBins); 
      return makeProjections( proj , argPack );
    }

    template <typename... ARGS> std::vector<TH1D*> makeProjections( const std::vector<Projection>& projections, const ARGS&... args )
    {
      return makeProjections( projections, ArgumentPack( args... ) );
    }

    template <typename... ARGS, 
              typename = std::enable_if_t< ! std::is_same<zeroType<ARGS...>, ArgumentPack>::value > > 
    TH1D* makeProjection( const Projection& projection, const ARGS&... args ) const
    {
      return makeProjection( projection, ArgumentPack(args...) );
    }

    template <typename... ARGS, 
              typename = std::enable_if_t< ! std::is_same<zeroType<ARGS...>, ArgumentPack>::value > > 
    TH2D* makeProjection( const Projection2D& projection, const ARGS&... args )
    {
      return makeProjection( projection, ArgumentPack(args...) );
    }

    template <typename functor> EventListSIMD& transform( functor&& fcn )
    {
      for ( auto& event : *this ) fcn( event );
      return *this;
    }
    static std::vector<float_v> makeEvent( const Event& event )
    {
      std::vector<float_v> rt( event.size() );
      for( unsigned i = 0 ; i != event.size(); ++i ) rt[i] = event[i];
      return rt;
    }
  };

} // namespace AmpGen
#endif

#endif
