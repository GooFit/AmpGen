#include "AmpGen/FastDT.h"
#include "AmpGen/BinDT.h"
#include <fstream>
#include <random>

using namespace AmpGen;

FastDT::FastDT( std::ifstream& stream, const unsigned& n_nodes )
{
  readFromStream(stream, n_nodes );
}

void FastDT::readFromStream( std::ifstream& stream, const int& n_nodes )
{
  std::string tmp;
  do 
  {
    std::getline(stream, tmp);
    if( tmp == "") break;
    auto tokens = split( tmp, ' ');
    if( tokens.size() != 4 ) FATAL("Could not read node: " << tmp );
    m_nodes.emplace_back( stoi( tokens[0]), stoi( tokens[1]), stoi(tokens[2]), stod( tokens[3] ) );
  } while( tmp != "");
}

void FastDT::serialise( std::ofstream& stream )
{
  for( auto& node : m_nodes ) stream << node.index << " " << node.left << " " << node.right << " " << node.cutValue << std::endl; 
}


FastDT::FastDT( const ArgumentPack& args )
{
  m_minEvents           = args.getArg<MinEvents>( 20 );
  m_maxDepth            = args.getArg<MaxDepth>( 999 );
  m_dim                 = args.getArg<Dim>( 5 );
  m_minStep             = std::vector<double>(m_dim,0);
}

int FastDT::findNode( const double* event ) const
{
  int address = m_nodes.size()-1;
  while( m_nodes[address].index != 0 )
  {
    auto node = m_nodes[address];
    if( address < 0 or address >= m_nodes.size() ) ERROR("Invalid node: " << address ); 
    // INFO("@ address: " << address << "next node = " << (event[node.index-1] > node.cutValue ? node.right : node.left)  );
   //  if( node.index -1 < 0 or node.index -1 >= m_dim ) ERROR("Invalid index: " << node.index -1 );
    address = event[node.index-1] > node.cutValue ? node.right : node.left; 
    if( address < 0 ) return std::abs(address) - 1; 
  }
  FATAL("Fallen out of DT ...");
  return m_nodes[address].left; 
}


double nnVariance(std::vector<double*>& evts, const unsigned& index)
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


void FastDT::refreshQueue(std::vector<double*>& evts, std::queue<unsigned>& indexQueue, const unsigned& depth)
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

std::pair<double,double> FastDT::bestCut_ls(const std::vector<double*>& source, 
                          const std::vector<double*>& target, 
                          int index, 
                          const size_t& dim,
                          const size_t& minEvents)
{
  auto w                       = [dim](const double* evt){ return *(evt+dim); };
  auto f                       = [index](const double* evt){ return *(evt+index); };
  double maxChi2               = -9999;
  double optPos                = 0;
  double s_w1      = 0;
  double s_wt      = parallel_accumulate(source.begin(), source.end()              , 0.0, w); //  sum_weights);
  double t_w1      = parallel_accumulate(target.begin(), target.begin() + minEvents, 0.0, w); //  sum_weights);
  double t_wt      = parallel_accumulate(target.begin(), target.end()  , 0.0, w);             //  sum_weights);
  unsigned int s_p1 = 0;
  unsigned int s_pt = source.size();
  auto sourceIt = source.begin();
  auto p0 = source[0][index];

  for(unsigned i = minEvents; i < target.size() - minEvents; i++ ) {
    t_w1  += w( target[i] );
    double pos = 0.5 * ( f(target[i]) + f(target[i + 1]) );
    for(; sourceIt != source.end(); ++sourceIt) {
      double positionOfThis = f(*sourceIt);
      if ( positionOfThis > pos ) break;
      s_w1  += w(*sourceIt);
      s_p1++;
    }
    double chi2 = (s_w1 - t_w1)*(s_w1 - t_w1)/(s_w1 + t_w1) + (s_wt - s_w1 - t_wt + t_w1)*(s_wt - s_w1 - t_wt + t_w1)/(s_wt - s_w1 + t_wt - t_w1);
    if ( chi2 > maxChi2 and s_p1 >= minEvents and (s_pt-s_p1) >= minEvents and pos - p0 > m_minStep[index] ) 
    { 
      maxChi2 = chi2;
      optPos  = pos;
    }
  }
  return {optPos, maxChi2};
}

int FastDT::makeNodes(std::vector<double*>& evts, std::queue<unsigned> indexQueue, const unsigned int& depth)
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
int FastDT::makeNodes( std::vector<double*> source, std::vector<double*> target)
{
  return makeNodes( source, target, makeQueue(), 0 );
}

std::vector<int> FastDT::makeQueue()
{
  DEBUG("Making queue ... " << m_dim );
  std::vector<int> iq( m_dim );
  std::iota( iq.begin(), iq.end() ,0 );
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(iq.begin(), iq.end(), g);
  return iq; 
}

int FastDT::makeNodes( std::vector<double*> source, std::vector<double*> target, std::vector<int> indexQueue, const unsigned int& depth)
{
  DEBUG("Reached depth: " << depth << " " << source.size() << " " << target.size() ); 
  if ( depth >= m_maxDepth || source.size() < 2 * m_minEvents || target.size() < 2 * m_minEvents ) {
    return --m_endNodeCounter;
  }
  unsigned bestIndex = 0; 
  double bestChi2 = 0; 
  double optPos = 0;
  if( std::all_of( indexQueue.begin(), indexQueue.end(), [](auto& arg){ return arg ==-1 ; } ) ) indexQueue = makeQueue();
  std::vector<double> chi2_v( m_dim ); 
  if( m_strategy == Strategy::best )
  { 
    auto it = indexQueue.begin();
    auto bestIt = it; 
    for(; it != indexQueue.end(); ++it )
    {
      auto i = *it;
      if( *it == -1 ) continue; 
      auto sorter = [&i]( auto& a, auto& b ) { return *( a + i ) < *( b + i ); };
      parallel_sort(source.begin(), source.end(), sorter, std::max(256ul, source.size() / 4) );
      parallel_sort(target.begin(), target.end(), sorter, std::max(256ul, target.size() / 4) );
      auto [op,chi2] = bestCut_ls(source, target, i, m_dim, m_minEvents);
      if( chi2 > bestChi2 ){
        bestChi2 = chi2; 
        bestIt = it;
        optPos = op; 
      } 
      chi2_v[i] = chi2; 
    }
    bestIndex = *it; 
    *it = -1;
    auto sorter = [bestIndex]( auto& a, auto& b ) { return *( a + bestIndex ) < *( b + bestIndex ); };
    parallel_sort(source.begin(), source.end(), sorter, std::max(256ul, source.size() / 4) );
    parallel_sort(target.begin(), target.end(), sorter, std::max(256ul, target.size() / 4) );
    DEBUG("Chi2_v = " << vectorToString( chi2_v, ", ") << " cutting @ [" << bestIndex << ", " << optPos ); 
  }
  else if( m_strategy == Strategy::random ) 
  {
    DEBUG("Selecting cut randomly from queue: " << vectorToString(indexQueue, ", " ) );
    auto it = std::find_if_not( std::begin(indexQueue), std::end(indexQueue), [](auto& it){ return it==-1;} );
    bestIndex = *it; 
    DEBUG("Cutting on: " << bestIndex );
    auto sorter = [bestIndex]( auto& a, auto& b ) { return *( a + bestIndex ) < *( b + bestIndex ); };
    parallel_sort(source.begin(), source.end(), sorter, std::max(256ul, source.size() / 4) );
    parallel_sort(target.begin(), target.end(), sorter, std::max(256ul, target.size() / 4) );
    auto [op,chi2] = bestCut_ls(source, target, bestIndex, m_dim, m_minEvents);
    bestChi2 = chi2; 
    optPos = op; 
    *it = -1;
  }
  if( bestChi2 == -9999 ) return --m_endNodeCounter; 
  auto sc = std::lower_bound( source.begin(), source.end(), optPos, [bestIndex](const auto& a, const auto& b){return b >= *(a+bestIndex); } );
  auto tc = std::lower_bound( target.begin(), target.end(), optPos, [bestIndex](const auto& a, const auto& b){return b >= *(a+bestIndex); } );
  std::vector<double*> sourceLeft(source.begin(), sc); 
  std::vector<double*> sourceRight(sc, source.end());
  std::vector<double*> targetLeft(target.begin(), tc);
  std::vector<double*> targetRight(tc, target.end());
  DEBUG("Cutting on " << bestIndex << " @ " << optPos << " chi2 = " << bestChi2 << " p[0] = " << *(source[0]+bestIndex)  << " s = [" << sourceLeft.size() << ", " << sourceRight.size() << "]; "
     << " t = [" << targetLeft.size() << ", " << targetRight.size() <<"]" );
  if ( sourceLeft.size() == 0 || sourceRight.size() == 0 ) {
    DEBUG( "No data to divide into block" );
    DEBUG( "Dividing:  source=(" << sourceLeft.size() << ", " << sourceRight.size() << "); target=("
        << targetLeft.size() << ", " << targetRight.size() << "); index=" << bestIndex
        << " cutting @ " << optPos << " zeroth source " << *(source[0] + bestIndex )  << " zeroth target: " << *(target[0]+bestIndex) 
        << " depth = " << depth );
    DEBUG( "Chi2 / index: " << vectorToString( chi2_v , ", ") );
    return --m_endNodeCounter;
  }
  auto left  = makeNodes( sourceLeft, targetLeft, indexQueue, depth + 1 );
  auto right = makeNodes( sourceRight, targetRight, indexQueue, depth + 1 );
  m_nodes.emplace_back( bestIndex+1, left, right, optPos );
  return m_nodes.size() -1;
}


