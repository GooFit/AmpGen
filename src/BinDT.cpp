#include "AmpGen/BinDT.h"

#include <string.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ratio>

#include "AmpGen/Utilities.h"
#include "AmpGen/Event.h"
#include <numeric>

using namespace AmpGen;

double nearestNeighbourVariance(std::vector<double*> evts, const size_t& index)
{
  auto dist  = [&index](const double* x, const double* y){ return fabs( *(x+index) -*(y+index)) ;};
  auto dist2 = [&index](const double* x, const double* y){ return ( *(x+index) -*(y+index))*( *(x+index) - *(y+index) ) ;};
  parallel_sort(evts.begin() 
                , evts.end() 
                , [index]( const double* evt1, const double* evt2 ) { return *( evt1 + index ) > *( evt2 + index ); }
                , std::max(256ul, evts.size() / 4) );
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

std::pair<std::vector<double*>, std::vector<double*>> splitEvents( std::vector<double*> evts, const size_t& index, double& midposition )
{
  auto sorter = [&index]( double* a, double* b ) { return *( a + index ) < *( b + index ); };
  std::sort( evts.begin(), evts.end(), sorter );
  unsigned int midpoint = evts.size() / 2;
  midposition =
    evts.size() % 2 == 0
    ? ( *( evts[midpoint - 1] + index ) + *( evts[midpoint] + index ) ) / 2
    : ( *( evts[midpoint + 1] + index ) + *( evts[midpoint - 1] + index ) + *( evts[midpoint] + index ) ) / 3;
  midpoint += midposition > *( evts[midpoint] + index );
  std::vector<double*> leftEvents( evts.begin(), evts.begin() + midpoint );
  std::vector<double*> rightEvents( evts.begin() + midpoint, evts.end() );
  return std::make_pair(leftEvents, rightEvents);
}

unsigned int BinDT::getBinNumber( const Event& evt ) const
{
  auto tEvent = m_functors( evt );
  auto leaf   = ( *m_top )( tEvent.data() );
  return leaf->voxNumber();
}

unsigned int BinDT::getBin( const Event& evt ) const
{
  auto tEvent = m_functors( evt );
  return ( *m_top )( tEvent.data() )->binNumber();
}

unsigned int BinDT::getBin( const double* evt ) const { return ( *m_top )( evt )->binNumber(); }

unsigned int BinDT::getBinNumber( const double* evt ) const { return ( *m_top )( evt )->voxNumber(); }

unsigned int BinDT::size() const { return m_endNodes.size(); }

std::function<std::vector<double>( const Event& )> BinDT::makeDefaultFunctors()
{
  if ( m_dim == 5 ) {
    return []( const Event& evt ) -> std::vector<double> {
      return {evt.s( 1, 2, 3 ), evt.s( 0, 1 ), evt.s( 0, 2 ), evt.s( 2, 3 ), evt.s( 0, 1, 2 )};
    };
  }
  if ( m_dim == 2 ) {
    DEBUG( "Problem has 2 d.o.f.s -> using Dalitz coordinates" );
    return []( const Event& evt ) -> std::vector<double> { return {evt.s( 0, 1 ), evt.s( 1, 2 )}; };
  }
  DEBUG( "No functors found for dim = " << m_dim );
  return nullptr;
}

BinDT::BinDT( const ArgumentPack& args )
{
  m_minEvents           = args.getArg<MinEvents>( 20 );
  m_maxDepth            = args.getArg<MaxDepth>( 999 );
  m_dim                 = args.getArg<Dim>( 5 );
  m_functors            = args.getArg<Functor>( makeDefaultFunctors() ).val;
  auto fname            = args.getArg<File>("").val;
  if ( fname != "" ) {
    auto stream = std::ifstream(fname);
    readFromStream(stream);
  }
}

void BinDT::readFromStream( std::istream& stream )
{
  std::map<std::string, std::pair<std::string, std::shared_ptr<INode>>> nodes;
  std::string line;
  std::string topAddress = "";
  while ( getline( stream, line ) ) {
    if( line == "end" ) break; 
    auto tokens = split( line, ' ' );
    try { 
      std::string address                = tokens[0];
      if ( topAddress == "" ) topAddress = address;
      if ( tokens.size() == 5 ) {
        int index      = stoi( tokens[1] );
        auto node      = std::make_shared<Decision>( index, stod( tokens[2] ), nullptr, nullptr );
        nodes[address] = std::make_pair( line, node );
      } else {
        auto node = std::make_shared<EndNode>( stoi( tokens[1] ) );
        if ( tokens.size() == 3 ) node->setBinNumber( stoi( tokens[2] ) );
        nodes[address] = std::make_pair( line, node );
        m_endNodes.push_back( node );
      }
    } catch(...){
      ERROR("Could not parse line: " << line );
    }
  }
  DEBUG( "Read: " << nodes.size()
     << " nodes from plain-text file" ); // , now re-establishing tree structure; EndNodes=" << edc );
  m_top = nodes[topAddress].second;
  for ( auto& node : nodes ) {
    Decision* node_ptr = dynamic_cast<Decision*>( node.second.second.get() );
    if ( node_ptr != nullptr ) {
      auto tokens = split( node.second.first, ' ' );
      node_ptr->setChildren( nodes[tokens[3]].second, nodes[tokens[4]].second );
    }
  }
}

void BinDT::serialize( std::ofstream& output )
{
  output << std::setprecision( 17 );
  m_top->serialize( output );
  output << "end";
}

void BinDT::serialize( const std::string& filename )
{
  std::ofstream output;
  output.open( filename );
  serialize( output );
  output.close();
}

std::shared_ptr<BinDT::INode> BinDT::makeNodes( std::vector<double*> evts )
{
  INFO( "Making nodes" );
  std::queue<unsigned int> iq;
  refreshQueue( evts, iq, 0 );
  auto tStartDT    = std::chrono::high_resolution_clock::now();
  auto node        = makeNodes( evts, iq, 0 );
  auto tNow        = std::chrono::high_resolution_clock::now();
  double timeTotal = std::chrono::duration<double, std::milli>( tNow - tStartDT ).count();
  INFO( "Time taken = " << timeTotal << " number of nodes = " << m_endNodes.size() );
  return node;
}

std::shared_ptr<BinDT::INode> BinDT::makeNodes( std::vector<double*> evts, std::queue<unsigned int> indexQueue,
    const unsigned int& depth )
{
  unsigned int index = indexQueue.front();
  DEBUG( "Depth = " << depth << " nodes = " << m_endNodes.size() << " maxDepth = " << m_maxDepth
      << " evts.size() = " << evts.size() << " minEvents = " << m_minEvents );
  if ( evts.size() < 2 * m_minEvents || depth > m_maxDepth ) {
    DEBUG( "Returning end node as " << evts.size() << " less than 2 x " << m_minEvents );
    auto node = std::make_shared<EndNode>( m_endNodes.size() );
    m_endNodes.push_back( node );
    return node;
  }
  double mid_point; 
  auto split_events = splitEvents(evts, index, mid_point);
  indexQueue.pop();
  if ( indexQueue.empty() ) refreshQueue( evts, indexQueue, depth + 1 );
  auto left  = makeNodes(split_events.first, indexQueue, depth + 1);
  auto right = makeNodes(split_events.second, indexQueue, depth + 1);
  auto node  = std::make_shared<Decision>( index, mid_point, right, left );
  return node;
}

std::shared_ptr<BinDT::INode> BinDT::makeNodes( std::vector<double*> source, std::vector<double*> target )
{
  std::queue<unsigned int> iq;
  refreshQueue( source, iq, 0 );
  m_top     = makeNodes(source, target, iq, 0);
  return m_top;
}
double bestCut_lineSearch(const std::vector<double*>& source, 
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
  double t_wt      = parallel_accumulate(target.begin(), target.end()  , 0.0, w); //  sum_weights);
  unsigned int s_p1 = 0;
  unsigned int s_pt = source.size();
  auto sourceIt = source.begin();
  for(size_t i = minEvents; i < target.size() - minEvents; i++ ) {
    t_w1  += w( target[i] );
    double pos = 0.5 * ( f(target[i]) + f(target[i + 1]) );
    for(; sourceIt != source.end(); ++sourceIt) {
      double positionOfThis = f(*sourceIt);
      if ( positionOfThis > pos ) break;
      s_w1  += w(*sourceIt);
      s_p1++;
    }
    //INFO( sourceLeftWeight << " " << targetLeftWeight << " " << sourceRightWeight << " " << targetRightWeight ); 
    double chi2 = (s_w1 - t_w1)*(s_w1 - t_w1)/(s_w1 + t_w1) + (s_wt - s_w1 - t_wt + t_w1)*(s_wt - s_w1 - t_wt + t_w1)/(s_wt - s_w1 + t_wt - t_w1);
    if ( chi2 > maxChi2 && s_p1 >= minEvents && (s_pt-s_p1) >= minEvents ) {
      maxChi2 = chi2;
      optPos  = pos;
    }
  }
  return optPos;
}
double BinDT::getBestPost( const std::vector<double*>& source, const std::vector<double*>& target, int index, bool verbose )
{
  double optChi2                     = -9999;
  double optPos                      = 0;
  double targetLeftWeight            = 0;
  double targetRightWeight           = target.size();
  double sourceLeftWeight            = 0;
  double sourceRightWeight           = 0;
  unsigned int sourceLeftPopulation  = 0;
  unsigned int sourceRightPopulation = source.size();
  for ( auto& event : source ) sourceRightWeight += *( event + m_dim );
  auto sourceIt = source.begin(); 
  INFO("target = " << targetLeftWeight << " " << targetRightWeight << " " << sourceLeftWeight << " " << sourceRightWeight );
  for ( unsigned int i = 0; i < target.size() - m_minEvents; i++ ) {
    targetLeftWeight += target[i][m_dim];
    targetRightWeight -= target[i][m_dim];
    if ( i < m_minEvents ) continue;
    double pos = ( target[i][index] + target[i + 1][index] ) / 2.;
    for ( ; sourceIt != source.end(); ++sourceIt ) {
      double positionOfThis = *( *sourceIt + index );
      if ( positionOfThis > pos ) break;
      sourceLeftWeight += *( *sourceIt + m_dim );
      sourceRightWeight -= *( *sourceIt + m_dim );
      sourceLeftPopulation++;
      sourceRightPopulation--;
    }
    double chi2 = ( sourceLeftWeight - targetLeftWeight ) * ( sourceLeftWeight - targetLeftWeight ) /
      ( sourceLeftWeight + targetLeftWeight ) +
      +( sourceRightWeight - targetRightWeight ) * ( sourceRightWeight - targetRightWeight ) /
      ( sourceRightWeight + targetRightWeight );
    if ( chi2 > optChi2 && sourceLeftPopulation >= m_minEvents && sourceRightPopulation >= m_minEvents ) {
      optChi2 = chi2;
      optPos  = pos;
    }
  }
  return optPos;
}

std::shared_ptr<BinDT::INode> BinDT::makeNodes( std::vector<double*> source, std::vector<double*> target,
    std::queue<unsigned int> indexQueue, const unsigned int& depth )
{
  unsigned int index = indexQueue.front();
  if ( source.size() == 0 || target.size() == 0 ) {
    ERROR( "No data to divide into block at opening!" );
  }
  if ( depth > m_maxDepth || source.size() < 2 * m_minEvents || target.size() < 2 * m_minEvents ) {
    auto node = std::make_shared<EndNode>( m_endNodes.size() );
    m_endNodes.push_back( node );
    return node;
  }
  auto sorter = [&index]( auto& a, auto& b ) { return *( a + index ) < *( b + index ); };
  parallel_sort(source.begin(), source.end(), sorter, std::max(256ul, source.size() / 4) );
  parallel_sort(target.begin(), target.end(), sorter, std::max(256ul, target.size() / 4) );
  double optPos = bestCut_lineSearch(source, target, index, m_dim, m_minEvents);
  auto sc = std::lower_bound( source.begin(), source.end(), optPos, [index](const auto& a, const auto& b){return b > *(a+index); } );
  auto tc = std::lower_bound( target.begin(), target.end(), optPos, [index](const auto& a, const auto& b){return b > *(a+index); } );
  std::vector<double*> sourceLeft(source.begin(), sc); 
  std::vector<double*> sourceRight(sc, source.end());
  std::vector<double*> targetLeft(target.begin(), tc);
  std::vector<double*> targetRight(tc, target.end());
  if ( sourceLeft.size() == 0 || sourceRight.size() == 0 ) {
    DEBUG( "No data to divide into block" );
    DEBUG( "Dividing:  source=(" << sourceLeft.size() << ", " << sourceRight.size() << "); target=("
        << targetLeft.size() << ", " << targetRight.size() << "); index=" << index
        << " depth = " << depth );
    auto node = std::make_shared<EndNode>( m_endNodes.size() );
    m_endNodes.push_back( node );
    return node;
  }
  indexQueue.pop();
  if ( indexQueue.empty() ) refreshQueue( source, indexQueue, depth + 1 );
  auto left  = makeNodes( sourceLeft, targetLeft, indexQueue, depth + 1 );
  auto right = makeNodes( sourceRight, targetRight, indexQueue, depth + 1 );
  auto node  = std::make_shared<Decision>( index, optPos, right, left );
  return node;
}

void BinDT::refreshQueue( const std::vector<double*>& evts, std::queue<unsigned int>& indexQueue,
    const unsigned int& depth )
{
  if ( evts.size() > m_minEvents * pow( 2, m_dim ) && depth + m_dim < m_maxDepth ) {
    if( m_queueOrdering.size() == 0 ) 
      for ( unsigned int i = 0; i < m_dim; ++i ) m_queueOrdering.push_back( i );
    for ( unsigned int i = 0; i < m_dim; ++i ) indexQueue.push( m_queueOrdering[i] );
  } else {
    std::vector<std::pair<unsigned int, double>> indices;
    for ( unsigned int i = 0; i < m_dim; ++i ) indices.emplace_back( i, nearestNeighbourVariance( evts, i ) );
    std::sort( indices.begin(), indices.end(),
        []( auto& it1, auto& it2 ) { return it1.second > it2.second;} );
    for ( auto& item : indices ) indexQueue.push( item.first );
  }
}

BinDT::Decision::Decision( const unsigned int& index, const double& value, std::shared_ptr<INode> left,
    std::shared_ptr<INode> right )
  : INode(), m_left( left ), m_right( right ), m_index( index ), m_value( value )
{
  if ( m_left != nullptr ) m_left->m_parent   = this;
  if ( m_right != nullptr ) m_right->m_parent = this;
}

void BinDT::Decision::serialize( std::ostream& stream ) const
{
  stream << this << " " << m_index << " " << m_value << " " << m_left.get() << " " << m_right.get() << std::endl;
  m_left->serialize( stream );
  m_right->serialize( stream );
}

void BinDT::Decision::setChildren( std::shared_ptr<INode> l, std::shared_ptr<INode> r )
{
  m_left            = l;
  m_right           = r;
  m_left->m_parent  = this;
  m_right->m_parent = this;
}

BinDT::EndNode::EndNode( const unsigned int& no, const unsigned int& binNumber ) : 
  m_voxNumber( no ), 
  m_binNumber( binNumber ) {}

const BinDT::EndNode* BinDT::EndNode::operator()( const double* evt ) const { return this; }

void BinDT::EndNode::serialize( std::ostream& stream ) const
{
  stream << this << " " << m_voxNumber << " " << m_binNumber << std::endl;
}

const BinDT::EndNode* BinDT::Decision::operator()( const double* evt ) const
{
  return *( evt + m_index ) < m_value ? ( *m_right )( evt ) : ( *m_left )( evt );
}

void BinDT::Decision::visit( const std::function<void(BinDT::INode*)>& visit_function )
{
  visit_function( this );
  m_left->visit( visit_function );
  m_right->visit( visit_function );
}
