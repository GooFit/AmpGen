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

#include "AmpGen/simd/iterator.h"
#include "AmpGen/simd/utils.h"
#include "AmpGen/Store.h"

namespace AmpGen
{
  class CompiledExpressionBase; 
  class EventListSIMD
  {
  private:
    Store<float_v, Alignment::AoS>   m_data      {};
    std::vector<float_v>             m_weights   {};
    std::vector<float_v>             m_genPDF    {};
    EventType                        m_eventType {};
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
    EventListSIMD( const EventList& other );     
    const float_v* data() const { return m_data.data(); }
    operator Store<float_v, Alignment::AoS> () const { return m_data ; }
    const auto& store()                     const { return m_data; }    
    const Event at(const unsigned& p)       const { return EventListSIMD::operator[](p) ; }
    const float_v* block(const unsigned& p) const { return m_data.data() + p * m_data.nFields(); }
          float_v* block(const unsigned& p)       { return m_data.data() + p * m_data.nFields(); }
    float_v weight(const unsigned& p) const { return m_weights[p]; }
    float_v genPDF(const unsigned& p) const { return m_genPDF[p]; }
    
    void setWeight( const unsigned& block, const float_v& w, const float_v& g=1)
    {
      m_weights[block] = w;
      m_genPDF[block] = g;
    } 
    void resize( const unsigned nEvents )
    {
      m_data = Store<float_v, Alignment::AoS>( nEvents, m_eventType.eventSize() );
      m_weights.resize( aligned_size(), 1);
      m_genPDF.resize( aligned_size(), 1 );
    }
    const Event operator[]( const size_t&) const;
    std::array<Event, utils::size<float_v>::value> scatter(unsigned) const;
    void gather(const std::array<Event, utils::size<float_v>::value>&, unsigned);   
    auto begin() const { return make_scatter_iterator<utils::size<float_v>::value>(0,this); }
    auto   end() const { return make_scatter_iterator<utils::size<float_v>::value>(size(), (const EventListSIMD*)(nullptr) ); } 
    auto begin()       { return make_scatter_iterator<utils::size<float_v>::value, true>(0, this); }
    auto   end()       { return make_scatter_iterator<utils::size<float_v>::value, true>(size(), (EventListSIMD*)(nullptr) ); }
    EventType eventType()                         const { return m_eventType; }
    size_t aligned_size()                         const { return m_data.aligned_size(); }
    double integral()                             const;
    size_t eventSize()                            const { return m_data.nFields(); }
    size_t size()                                 const { return m_data.size(); }
    size_t nBlocks()                              const { return m_data.nBlocks(); }
    void setEventType( const EventType& type ) { m_eventType = type; }
    void add( const EventListSIMD& evts );
    void loadFromTree( TTree* tree, const ArgumentPack& args ); 
    void loadFromFile( const std::string& fname, const ArgumentPack& args );
    void clear();

    TTree* tree( const std::string& name, const std::vector<std::string>& extraBranches = {} ) const;
    

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
