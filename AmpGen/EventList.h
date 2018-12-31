#ifndef AMPGEN_EVENTLIST_H
#define AMPGEN_EVENTLIST_H

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"

#include <chrono>
#include <functional>
#include <numeric>
#include <cstddef>
#include <algorithm>

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#ifdef __USE_OPENMP__
#include <omp.h>
#endif

namespace AmpGen
{
  DECLARE_ARGUMENT_DEFAULT( Bins, size_t, 100 );

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
    EventList() = default;
    EventList( const EventType& type );
    template < class ... ARGS > EventList( const std::string& fname, const EventType& evtType, const ARGS&... args ) : EventList(evtType) 
    {
      loadFromFile( fname, ArgumentPack(args...) );
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
    size_t registerExpression( const CompiledExpressionBase& expression, const size_t& size_of =0 );

    template <class FUNCTOR>
    unsigned int extendEvent( const std::string& name, FUNCTOR func )
    {
      unsigned int index = this->begin()->size();
      for ( auto& evt : *this ) {
        evt.extendEvent( func( evt ) );
      }
      m_extensions[name] = index;
      return index;
    }
    template <class FCN>
    void updateCache( const FCN& fcn, const size_t& index )
    {
      #pragma omp parallel for
      for ( unsigned int i = 0; i < size(); ++i ) {
        ( *this )[i].setCache( fcn( getEvent(i) ), index );
      }
    }
    TH2D* makeProjection( const Projection2D& projection, const ArgumentPack& args );

    std::vector<TH1D*> makePlots( const std::vector<Projection>& projections, const ArgumentPack& args );

    template <class... ARGS>
    std::vector<TH1D*> makeDefaultPlots( const ARGS&... args )
    {
      auto argPack = ArgumentPack( args... );
      size_t nBins = argPack.getArg<Bins>(100);
      return makePlots( eventType().defaultProjections(nBins), argPack );
    }
    template <class... ARGS> 
    std::vector<TH1D*> makePlots( const std::vector<Projection>& projections, const ARGS&... args )
    {
      return makePlots( projections, ArgumentPack( args... ) );
    }
    template <class... ARGS>
    TH1D* makeProjection( const Projection& projection, const ARGS&... args )
    {
      return makeProjection( projection, ArgumentPack(args...) );
    }
    TH1D* makeProjection( const Projection& projection, const ArgumentPack& args );
    
    template <class... ARGS>
    TH2D* makeProjection( const Projection2D& projection, const ARGS&... args )
    {
      return makeProjection( projection, ArgumentPack(args...) );
    }
    template <class FCN>
    EventList& transform( FCN&& fcn )
    {
      for ( auto& event : m_data ) fcn( event );
      return *this;
    }
    template <class FCN>
    void filter( FCN&& fcn ){
      size_t currentSize = size();
      m_data.erase( std::remove_if( m_data.begin(), m_data.end(), fcn ) , m_data.end() );
      INFO("Filter removes: " << currentSize - size() << " / " << currentSize << " events");
    }
  };
  DECLARE_ARGUMENT( Prefix, std::string );
  DECLARE_ARGUMENT( LineColor, int );
  DECLARE_ARGUMENT( DrawStyle, std::string );
  DECLARE_ARGUMENT( Selection, std::function<bool( const Event& )> );
  DECLARE_ARGUMENT( WeightFunction, std::function<double( const Event& ) > );
  DECLARE_ARGUMENT( Branches, std::vector<std::string> );
  DECLARE_ARGUMENT( EntryList, std::vector<size_t> );
  DECLARE_ARGUMENT_DEFAULT( GetGenPdf, bool, false );
  DECLARE_ARGUMENT_DEFAULT( CacheSize, size_t , 0 );
  DECLARE_ARGUMENT_DEFAULT( Filter, std::string , "");
  DECLARE_ARGUMENT_DEFAULT( WeightBranch, std::string, "" );      
  DECLARE_ARGUMENT_DEFAULT( ApplySym, bool, 0 );  
} // namespace AmpGen

#endif
