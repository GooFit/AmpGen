#ifndef AMPGEN_FASTDT_H
#define AMPGEN_FASTDT_H
#include <memory.h>
#include <stddef.h>
#include <array>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
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

namespace AmpGen {
  class FastDT
  {
    public:
      enum Strategy { best, random }; 
      struct Node {
        int index; 
        int left; 
        int right; 
        double   cutValue; 
        Node() = default;
        Node(const int& index,
            const int& left, 
            const int& right,
            const double& cutValue) : index(index), left(left), right(right), cutValue(cutValue) {};
      };
      FastDT() = default; 
      template <typename ...arg_types> FastDT( const arg_types&... args) 
      : FastDT( ArgumentPack(args...) ) {}
      FastDT( const ArgumentPack& args );
      FastDT( const std::string& textFile );
      FastDT( std::ifstream&, const unsigned&);
      int getBinNumber( const double* event) const { return findNode(event) ; }
      int findNode(const double* event) const;
      std::pair<double,double> bestCut_ls(const std::vector<double*>& source, 
                          const std::vector<double*>& target, 
                          int index, 
                          const size_t& dim,
                          const size_t& minEvents); 
      std::vector<Node> const_nodes() const { return m_nodes; }
      std::vector<Node> m_nodes;
      std::queue<unsigned> m_queueOrdering; 
      unsigned     m_dim            = {0};
      unsigned     m_minEvents      = {0};
      unsigned     m_maxDepth       = {0};
      int          m_endNodeCounter = {0};       
      Strategy     m_strategy       = {Strategy::random};
      std::vector<double> m_minStep; 
      std::vector<int> makeQueue(); 
      void setStep( const std::vector<double>& step ) { m_minStep = step; }
      int makeNodes(std::vector<double*>&, std::queue<unsigned>, const unsigned&);
      int makeNodes(std::vector<double*>, std::vector<double*>);
      int makeNodes(std::vector<double*>, std::vector<double*>, std::vector<int>, const unsigned&);
      void refreshQueue(std::vector<double*>& evts, std::queue<unsigned>& indexQueue, const unsigned& depth);
      void serialise( std::ofstream& );
      void setQueueOrdering( std::vector<unsigned>& ){};
      void readFromStream( std::ifstream& stream, const int& n_nodes );
  };
}
#endif
