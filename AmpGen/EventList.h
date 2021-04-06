#ifndef AMPGEN_EVENTLIST_H
#define AMPGEN_EVENTLIST_H

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/Units.h"

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
  namespace PlotOptions { DECLARE_ARGUMENT(Bins, size_t); }
  class CompiledExpressionBase; 
  class EventList
  {
  private:
    std::vector<Event>                  m_data              = {};
    EventType                           m_eventType         = {};
    std::map<std::string, unsigned int> m_extensions        = {};
        
  public:
    typedef Event value_type;
    EventList() = default;
    EventList( const EventType& type );
    template < class ... ARGS > EventList( const std::string& fname, const EventType& evtType, const ARGS&... args ) : EventList(evtType) 
    {
      loadFromFile( fname, ArgumentPack(args...) );
    }
    template < class ... ARGS > EventList( const std::string& fname, const ARGS&... args ) : EventList() 
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
    const EventList& store()                      const { return *this;}    
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
    size_t aligned_size()                         const { return m_data.size() ; }
    size_t nBlocks()                              const { return m_data.size() ; }
    double integral()                             const;
    const double* block(const unsigned pos) const { return m_data[pos].address(); }
    real_t weight( const size_t& pos)             const { return m_data[pos].weight(); }
    real_t genPDF( const size_t& pos)             const { return m_data[pos].genPdf(); }
    unsigned key(const std::string& key)          const 
    { 
      auto it = m_extensions.find(key);
      if( it == m_extensions.end() ) return m_data[0].size() - 1;
      return it->second; 
    }
    void reserve( const size_t& size );
    void resize ( const size_t& size );
    void push_back( const Event& evt );
    void emplace_back( const Event& evt);
    void setEventType( const EventType& type ) { m_eventType = type; }
    void add( const EventList& evts );
    void loadFromTree( TTree* tree, const ArgumentPack& args ); 
    void loadFromFile( const std::string& fname, const ArgumentPack& args );
    void clear();
    void setWeight( const unsigned int& pos, const double& w, const double&g=+1)
    {
      m_data[pos].setWeight(w);
      m_data[pos].setGenPdf(g);
    }
    void setGenPDF( const unsigned int& pos, const double& g)
    {
      m_data[pos].setGenPdf(g);
    }
    void erase( const std::vector<Event>::iterator& begin, const std::vector<Event>::iterator& end );

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

    template <typename functor> EventList& transform( functor&& fcn )
    {
      for ( auto& event : m_data ) fcn( event );
      return *this;
    }
    
    template <typename functor> void filter( functor&& fcn )
    {
      unsigned currentSize = size();
      m_data.erase( std::remove_if( m_data.begin(), m_data.end(), fcn ) , m_data.end() );
      INFO("Filter retains " << size() << " / " << currentSize << " events");
    }

    template <typename functor> unsigned count( functor&& fcn ) const 
    {
      return std::count_if( std::begin(*this), std::end(*this), fcn );
    }
  }; 
  DECLARE_ARGUMENT(Branches, std::vector<std::string>);      /// Branch names containing kinematic information 
  DECLARE_ARGUMENT(ExtraBranches, std::vector<std::string>); /// additional information about the event to include
  DECLARE_ARGUMENT(IdBranches, std::vector<std::string>);    /// Branches containing PID information, used if the names of particles are incorrect (looking at you, DTF)  
  DECLARE_ARGUMENT(EntryList, std::vector<size_t>);
  DECLARE_ARGUMENT(GetGenPdf, bool);
  DECLARE_ARGUMENT(Filter, std::string);
  DECLARE_ARGUMENT(WeightBranch, std::string);      
  DECLARE_ARGUMENT(ApplySym, bool);  
  DECLARE_ARGUMENT(WeightFunction, std::function<double( const Event& )>);
  DECLARE_ARGUMENT(InputUnits, AmpGen::Units);
} // namespace AmpGen

#endif
