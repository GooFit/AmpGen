#ifndef AMPGEN_EVENTLIST_H
#define AMPGEN_EVENTLIST_H

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"

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

namespace AmpGen
{

  DECLARE_ARGUMENT(Bins, size_t);

  class CompiledExpressionBase; 
  class EventList
  {
  private:
    std::vector<Event>                  m_data              = {};
    EventType                           m_eventType         = {};
    std::map<uint64_t, unsigned int>    m_pdfIndex          = {};
    std::map<std::string, unsigned int> m_extensions        = {};
    double                              m_norm              = {0};
    size_t                              m_lastCachePosition = {0}; 
  public:
    typedef Event value_type;
    EventList() = default;
    EventList( const EventType& type );
    template < class ... ARGS > EventList( const std::string& fname, const EventType& evtType, const ARGS&... args ) : EventList(evtType) 
    {
      loadFromFile( fname, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventList( const std::vector<std::string>& fname, const EventType& evtType, const ARGS&... args ) : EventList(evtType) 
    {
      for( auto& f : fname ) loadFromFile( f, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventList( TTree* tree, const EventType& evtType, const ARGS&... args ) : EventList(evtType)
    {
      loadFromTree( tree, ArgumentPack(args...) );
    }
    
    void resetCache();
    std::vector<Event>::reverse_iterator rbegin()       { return m_data.rbegin(); }
    std::vector<Event>::reverse_iterator rend()         { return m_data.rend(); }
    std::vector<Event>::iterator begin()                { return m_data.begin(); }
    std::vector<Event>::iterator end()                  { return m_data.end(); }
    Event& operator[]( const size_t& pos )              { return m_data[pos]; }
    real_t* getEvent( const size_t& index )             { return (( *this )[index] ); }
    const real_t* getEvent( const size_t& index ) const { return (const real_t*)(( *this )[index] ); }
    std::vector<Event>::const_iterator begin()    const { return m_data.cbegin(); }
    std::vector<Event>::const_iterator end()      const { return m_data.cend(); }
    const Event& operator[]( const size_t& pos )  const { return m_data[pos]; }
    EventType eventType()                         const { return m_eventType; }
    const Event& at( const size_t& pos )          const { return m_data[pos]; }
    size_t size()                                 const { return m_data.size(); }
    double integral()                             const;
    double norm();

    void reserve( const size_t& size ) { m_data.reserve( size ); }
    void push_back( const Event& evt ) { m_data.push_back( evt ); }
    void setEventType( const EventType& type ) { m_eventType = type; }
    void add( const EventList& evts );
    void loadFromTree( TTree* tree, const ArgumentPack& args ); 
    void loadFromFile( const std::string& fname, const ArgumentPack& args );
    void printCacheInfo( const unsigned int& nEvt = 0 );
    void clear();
    void erase( const std::vector<Event>::iterator& begin, const std::vector<Event>::iterator& end );

    TTree* tree( const std::string& name, const std::vector<std::string>& extraBranches = {} );
    
    size_t getCacheIndex( const CompiledExpressionBase& PDF, bool& status ) const;
    size_t getCacheIndex( const CompiledExpressionBase& PDF ) const;
    template <class T>
    size_t registerExpression(const T& expression, const size_t& size_of=0)
    {
      auto key = FNV1a_hash( expression.name() );
      auto pdfIndex = m_pdfIndex.find( key );
      if ( pdfIndex != m_pdfIndex.end() ) {
        return pdfIndex->second;
      } else {
        size_t lcp            = m_lastCachePosition;
        size_t expression_size = size_of == 0 ? 
          expression.returnTypeSize() / sizeof(complex_t) : size_of; 
        if (lcp >= at( 0 ).cacheSize() ) { 
          WARNING("Cache index " << lcp << " exceeds cache size = " 
                                 << at(0).cacheSize() << " resizing to " 
                                 << lcp + expression_size );
          resizeCache( lcp + expression_size );
        }
        m_pdfIndex[key] = m_lastCachePosition;
        m_lastCachePosition += expression_size; 
        return lcp;
      }
    }

    template <class FCN>
    void updateCache( const FCN& fcn, const size_t& index )
    {
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for ( unsigned int i = 0; i < size(); ++i ) {
        ( *this )[i].setCache(fcn(getEvent(i)), index);
      }
    }
    void reserveCache(const size_t& index);
    void resizeCache(const size_t& newCacheSize );
    TH1D* makeProjection(const Projection& projection  , const ArgumentPack& args = ArgumentPack()) const; 
    TH2D* makeProjection(const Projection2D& projection, const ArgumentPack& args = ArgumentPack()) const;
    std::vector<TH1D*> makeProjections( const std::vector<Projection>& projections, const ArgumentPack& args );

    template <class... ARGS> std::vector<TH1D*> makeDefaultProjections( const ARGS&... args )
    {
      auto argPack = ArgumentPack( args... );
      size_t nBins = argPack.getArg<Bins>(100);
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

    template <class FCN> EventList& transform( FCN&& fcn )
    {
      for ( auto& event : m_data ) fcn( event );
      return *this;
    }
    
    template <class FCN> void filter( FCN&& fcn ){
      size_t currentSize = size();
      m_data.erase( std::remove_if( m_data.begin(), m_data.end(), fcn ) , m_data.end() );
      INFO("Filter removes: " << currentSize - size() << " / " << currentSize << " events");
    }
  };
  DECLARE_ARGUMENT(LineColor, int);
  DECLARE_ARGUMENT(DrawStyle, std::string);
  DECLARE_ARGUMENT(Selection, std::function<bool( const Event& )>);
  DECLARE_ARGUMENT(WeightFunction, std::function<double( const Event& )>);
  DECLARE_ARGUMENT(Branches, std::vector<std::string>);
  DECLARE_ARGUMENT(EntryList, std::vector<size_t>);
  DECLARE_ARGUMENT(GetGenPdf, bool);
  DECLARE_ARGUMENT(CacheSize, size_t);
  DECLARE_ARGUMENT(Filter, std::string);
  DECLARE_ARGUMENT(WeightBranch, std::string);      
  DECLARE_ARGUMENT(ApplySym, bool);  
  DECLARE_ARGUMENT(Prefix, std::string);
} // namespace AmpGen

#endif
