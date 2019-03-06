#ifndef AMPGEN_FASTDT_H
#define AMPGEN_FASTDT_H
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

namespace AmpGen {
  class FastDT
  {
    public:
      struct Node {
        unsigned int index; 
        int left; 
        int right; 
        double       cutValue; 
        Node(){};
        Node(const unsigned int& index,
            const int& left, 
            const int& right,
            const double& cutValue) : index(index), left(left), right(right), cutValue(cutValue) {};
      };
      FastDT() = default; 
      FastDT( const std::string& textFile );
      int findNode(const double* event);

      std::vector<Node> m_nodes;
      std::function<std::vector<double>(const Event&)> m_functors;
      std::queue<size_t> m_queueOrdering; 
      size_t       m_dim;
      size_t       m_minEvents;
      size_t       m_maxDepth       = {0};
      int          m_endNodeCounter = {0};       
      int makeNodes(std::vector<double*>& evts, std::queue<size_t> indexQueue, const size_t& depth);
      void refreshQueue(std::vector<double*>& evts, std::queue<size_t>& indexQueue, const size_t& depth);
  };
}
#endif
