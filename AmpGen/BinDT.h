#ifndef AMPGEN_BINDT_H
#define AMPGEN_BINDT_H
#include <memory.h>
#include <stddef.h>
#include <array>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class Event;

  DECLARE_ARGUMENT(MaxDepth, size_t );
  DECLARE_ARGUMENT(MinEvents, size_t );
  DECLARE_ARGUMENT(Dim, size_t );
  DECLARE_ARGUMENT(Functor, std::function<std::vector<real_t>( const Event& )>);
  DECLARE_ARGUMENT(File, std::string);

  class BinDT
  {
    public:
      class EndNode;

      class INode
      {
        public:
          INode() = default;
          virtual ~INode() = default;
          virtual const EndNode* operator()( const double* evt ) const = 0;
          virtual void serialize( std::ostream& stream ) const         = 0;
          virtual void visit( const std::function<void( INode* )>& visit_function ) = 0; 
          INode* m_parent = {nullptr};
      };
      class EndNode : public INode
    {
      public:
        EndNode( const unsigned int& no, const unsigned int& binNumber = 999 );
        const EndNode* operator()( const double* evt ) const override ;
        void serialize( std::ostream& stream ) const override;
        unsigned int voxNumber() const { return m_voxNumber; }
        unsigned int binNumber() const { return m_binNumber; }
        void setBinNumber( const unsigned int& binNumber ) { m_binNumber = binNumber; }
        void visit( const std::function<void(INode*)>& visit_function ) override { visit_function( this ); }
        friend class BinDT;
      private:
        unsigned int m_voxNumber;
        unsigned int m_binNumber;
    };

    class Decision : public INode
    {
      public:
        Decision( const unsigned int& index, const double& value, std::shared_ptr<INode> left = nullptr,
            std::shared_ptr<INode> right = nullptr );
        const EndNode* operator()( const double* evt ) const override;
        void serialize( std::ostream& stream ) const override;
        void setChildren( std::shared_ptr<INode> l, std::shared_ptr<INode> r );
        void visit( const std::function<void(INode*)>& visit_function ) override;
        friend class BinDT;
      private :
        std::shared_ptr<INode> m_left;
        std::shared_ptr<INode> m_right;
        unsigned int m_index;
        double m_value;
    };

    public:
      template <class... ARGS>
        BinDT( const ARGS&... args ) : BinDT( ArgumentPack( args... ) )
      {
      }
      template <class... ARGS>
        BinDT( const std::vector<double*>& addr, const ARGS&... args ) : BinDT( ArgumentPack( args... ) )
      {
        m_top = makeNodes( addr );
      }
      template <class... ARGS>
        BinDT( const EventList& events, const ARGS&... args ) : BinDT( ArgumentPack( args... ) )
      {
        std::vector<double> data( m_dim * events.size() );
        std::vector<double*> addresses( events.size() );
        size_t counter = 0;
        for ( auto& evt : events ) {
          auto val = m_functors( evt );
          for ( unsigned int i = 0; i < m_dim; ++i ) data[m_dim * counter + i] = val[i];
          addresses[counter]                                                   = &( data[m_dim * counter] );
          counter++;
        }
        INFO( "Making nodes" );
        m_top = makeNodes( addresses );
      }
      BinDT( const ArgumentPack& args );
      BinDT() = default;

      std::shared_ptr<INode> top() { return m_top; }
      double nnUniformity( std::vector<double*> evts, const unsigned int& index ) const;
      unsigned int getBinNumber( const Event& evt ) const;
      unsigned int getBinNumber( const double* evt ) const;
      unsigned int getBin( const Event& evt ) const;
      unsigned int getBin( const double* evt ) const;
      unsigned int size() const;
      void readFromStream( std::istream& stream );
      void serialize( std::ofstream& output );
      void serialize( const std::string& filename );
      void setQueueOrdering( const std::vector<size_t>& queueOrdering ){ m_queueOrdering = queueOrdering ; }
      std::vector<std::shared_ptr<EndNode>>& nodes() { return m_endNodes; }
      const std::vector<std::shared_ptr<EndNode>>& const_nodes() const { return m_endNodes; }
      std::vector<std::shared_ptr<EndNode>>::iterator begin() { return m_endNodes.begin(); }
      std::vector<std::shared_ptr<EndNode>>::iterator end() { return m_endNodes.end(); }

      std::function<std::vector<double>( const Event& )> makeDefaultFunctors();
      void refreshQueue( const std::vector<double*>& evts, std::queue<unsigned int>& indexQueue,
          const unsigned int& depth );
      std::shared_ptr<INode> makeNodes( std::vector<double*> evts, std::queue<unsigned int> indexQueue,
          const unsigned int& depth );
      std::shared_ptr<INode> makeNodes( std::vector<double*> evts );
      std::shared_ptr<INode> makeNodes( std::vector<double*> source, std::vector<double*> target );

      std::shared_ptr<INode> makeNodes( std::vector<double*> source, std::vector<double*> target,
          std::queue<unsigned int> indexQueue, const unsigned int& depth );
      void setFunctor( const std::function<std::vector<double>( const Event& )>& functors ) { m_functors = functors; }
    
    private:
      std::shared_ptr<INode> m_top;
      unsigned int m_dim;
      std::vector<std::shared_ptr<EndNode>> m_endNodes;
      std::function<std::vector<double>( const Event& )> m_functors;
      unsigned int m_minEvents;
      unsigned int m_maxDepth;
      std::vector<size_t> m_queueOrdering; 
      double getBestPost(const std::vector<double*>& source, const std::vector<double*>& target, int index, bool verbose = false );
  };

} // namespace AmpGen
#endif
