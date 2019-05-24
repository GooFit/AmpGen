#include "AmpGen/FastDT.h"
#include "AmpGen/BinDT.h"
#include <fstream>

using namespace AmpGen;

int FastDT::findNode( const double* event )
{
  int address = m_nodes.size()-1;
  while( m_nodes[address].index != 0 )
  {
    auto node = m_nodes[address];
    address = event[node.index-1] > node.cutValue ? node.right : node.left; 
    if( address < 0 ) return std::abs(address) - 1; 
  }
  return m_nodes[address].left; 
}

FastDT::FastDT( const std::string& fname ){
  //std::ifstream s(fname);
}

double nnVariance(std::vector<double*>& evts, const size_t& index)
{
  auto dist  = [&index](const double* x, const double* y){ return fabs( *(x+index) -*(y+index)) ;};
  auto dist2 = [&index](const double* x, const double* y){ return ( *(x+index) -*(y+index))*( *(x+index) - *(y+index) ) ;};
  std::sort( evts.begin(), evts.end(), [index]( const double* evt1, const double* evt2 ) { return *( evt1 + index ) > *( evt2 + index ); } );
  size_t size                 = evts.size();
  double x   = dist(evts[0], evts[1]) + dist(evts[size-1], evts[size-2]); 
  double x2  = dist2(evts[0], evts[1]) + dist2(evts[size-1], evts[size-2]); 
  for ( size_t i = 1; i < size-1; ++i ) {
    double low  = dist(evts[i], evts[i-1]); 
    double high = dist(evts[i], evts[i+1]);
    x  += low > high ? high : low;
    x2 += low > high ? high*high : low*low;
  }
  return x2/double(size) - x*x / double(size*size);
}

std::pair<std::vector<double*>, std::vector<double*>> _splitEvents(std::vector<double*>& evts, const size_t& index, double& midposition)
{
  auto sorter = [&index]( double* a, double* b ) { return *( a + index ) < *( b + index ); };
  std::sort( evts.begin(), evts.end(), sorter );
  unsigned int midpoint = evts.size() / 2;
  midposition =
    evts.size() % 2 == 0
    ? ( *( evts[midpoint-1] + index ) + *( evts[midpoint] + index ) ) / 2
    : ( *( evts[midpoint+1] + index ) + *( evts[midpoint-1] + index ) + *( evts[midpoint] + index ) ) / 3;
  midpoint += midposition > *( evts[midpoint] + index );
  std::vector<double*> leftEvents( evts.begin(), evts.begin() + midpoint );
  std::vector<double*> rightEvents( evts.begin() + midpoint, evts.end() );
  return std::make_pair(leftEvents, rightEvents);
}


void FastDT::refreshQueue(std::vector<double*>& evts, std::queue<size_t>& indexQueue, const size_t& depth)
{
  if ( evts.size() > m_minEvents * pow(2, m_dim) && depth + m_dim < m_maxDepth ) {
    if( m_queueOrdering.size() == 0 ){ for ( size_t i = 0; i < m_dim; ++i ) m_queueOrdering.push(i); }
    indexQueue = m_queueOrdering; 
  } 
  else 
  {
    std::vector<std::pair<unsigned int, double>> indices;
    for (size_t i = 0; i < m_dim; ++i) indices.emplace_back( i, nnVariance( evts, i ) );
    std::sort( indices.begin(), indices.end(), []( auto& it1, auto& it2 ) { return it1.second > it2.second;} );
    for ( auto& item : indices ) indexQueue.push( item.first );
  }
}

int FastDT::makeNodes(std::vector<double*>& evts, std::queue<size_t> indexQueue, const size_t& depth)
{
  if( depth == 0 && indexQueue.empty() ) refreshQueue(evts, indexQueue, depth + 1);
  size_t index = indexQueue.front();
  if ( evts.size() < 2 * m_minEvents || depth > m_maxDepth ){
    return --m_endNodeCounter; 
  }
  double cutValue; 
  auto split_events = _splitEvents(evts, index, cutValue);
  indexQueue.pop();
  if( indexQueue.empty() ) refreshQueue(evts, indexQueue, depth + 1);
  auto left  = makeNodes(split_events.first , indexQueue, depth + 1);
  auto right = makeNodes(split_events.second, indexQueue, depth + 1);
  m_nodes.emplace_back(index+1, left, right, cutValue);
  return m_nodes.size() -1;
}
