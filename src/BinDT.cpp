#include "AmpGen/BinDT.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <string.h>
#include <type_traits>

#include "AmpGen/Utilities.h"

using namespace AmpGen;

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

namespace AmpGen
{
  struct PackedDecision {
    uint8_t header;   //      = 2 bytes
    uint8_t index;    //      = 8 bytes
    uint32_t address; //      = 2 + 4 bytes
    uint32_t left;    //      = 18
    uint32_t right;   //      = 24 bytes
    double cutVal;    //      = 14 bytes
  };
  struct PackedNode {
    uint8_t header;
    uint8_t packingBit[1];
    uint32_t address;
    uint32_t binNumber;
    uint32_t voxNumber;
    uint16_t packingBits[3];
  };
} // namespace AmpGen
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
  std::ifstream* stream = args.getArg<Stream>().val;
  auto fname            = args.getArg<File>();
  if ( stream != nullptr && stream->is_open() ) {
    // if( readType & std::ios::binary )
    readFromBinary( *stream );
    // else readFromStream( *stream );
  } else if ( fname.name != "" ) {
    std::ifstream stream( fname.name, fname.mode | std::ios::in );
    if ( fname.mode & std::ios::binary )
      readFromBinary( stream );
    else
      readFromStream( stream );
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
  INFO( "Read: " << nodes.size()
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

void BinDT::readFromBinary( std::ifstream& stream )
{
  std::map<uint32_t, std::shared_ptr<INode>> nodes;
  std::map<uint32_t, std::pair<uint32_t, uint32_t>> fti;
  char buffer[PACKET_SIZE]; // 24 byte word for the nodes //
  size_t edc = 0;
  while ( stream.read( buffer, PACKET_SIZE ) ) {
    if ( ( uint8_t( buffer[0] ) ) == 0xAA ) {
      PackedDecision pd;
      memcpy( &pd, buffer, PACKET_SIZE );
      nodes.emplace( pd.address, std::make_shared<Decision>( pd.index, pd.cutVal ) );
      DEBUG( "header = " << std::hex << (uint8_t)pd.header << std::dec << " index = " << pd.index
          << " cut-val = " << pd.cutVal << " left-node = " << pd.left << " right-node=" << std::hex
          << pd.right << " " << std::dec << sizeof( PackedDecision ) << " " << PACKET_SIZE );
      fti[pd.address] = std::make_pair( pd.left, pd.right );
    } else if ( uint8_t( buffer[0] ) == 0xBB ) {
      edc++;
      PackedNode pd;
      memcpy( &pd, buffer, PACKET_SIZE );
      nodes.emplace( pd.address, std::make_shared<EndNode>( pd.voxNumber, pd.binNumber ) );
    } else {
      ERROR( "Packet header: " << std::hex << uint8_t( buffer[0] ) << " not recognised" );
    }
    if ( nodes.size() == 1 ) m_top = nodes.begin()->second;
  }
  INFO( "Read: " << nodes.size() << " nodes from binary file, now re-establishing tree structure; EndNodes=" << edc );

  for ( auto& node : nodes ) {
    Decision* node_ptr = dynamic_cast<Decision*>( node.second.get() );
    if ( node_ptr == nullptr ) continue;
    auto child_nodes = fti.find( node.first );
    if ( child_nodes == fti.end() ) {
      ERROR( "Child nodes not found!" );
    }
    auto pL = nodes.find( child_nodes->second.first );
    auto pR = nodes.find( child_nodes->second.second );
    if ( pL == nodes.end() || pR == nodes.end() ) {
      ERROR( "Nodes for " << node_ptr << " not found! success counter "
          << " " << child_nodes->second.first << " " << child_nodes->second.second );
    } else {
      node_ptr->setChildren( pL->second, pR->second );
    }
  }
}

void BinDT::writeToBinary( std::ofstream& stream )
{
  AddressCompressor counter;
  m_top->visit( [&stream, &counter]( const INode* node ) {
      std::array<char, PACKET_SIZE> word;
      const BinDT::Decision* decision = dynamic_cast<const BinDT::Decision*>( node );
      if ( decision != nullptr ) {
      PackedDecision pd;
      pd.header  = 0xAA;
      pd.address = counter[decision]; // uint32_t ( 0xFFFFFFFF &  (uint64_t)this );
      pd.index   = decision->m_index;
      pd.cutVal  = decision->m_value;
      pd.left    = counter[decision->m_left.get()];  // (uint32_t)m_left.get();
      pd.right   = counter[decision->m_right.get()]; // (uint32_t)m_right.get();
      memcpy( word.data(), &pd, PACKET_SIZE );
      } else if ( dynamic_cast<const BinDT::EndNode*>( node ) != nullptr ) {
      const BinDT::EndNode* endNode = dynamic_cast<const BinDT::EndNode*>( node );
      PackedNode pd;
      pd.header    = 0xBB;
      pd.address   = counter[endNode];
      pd.binNumber = endNode->m_binNumber;
      pd.voxNumber = endNode->m_voxNumber;
      memcpy( word.data(), &pd, PACKET_SIZE );
      }
      stream.write( (char*)( &word ), sizeof( word ) );
  } );
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

double BinDT::nnUniformity( std::vector<double*> evts, const unsigned int& index ) const
{

  std::sort( evts.begin(), evts.end(),
      [index]( const double* evt1, const double* evt2 ) { return *( evt1 + index ) > *( evt2 + index ); } );
  double averageNNdistance = 0;
  int size                 = evts.size();
  for ( int i = 0; i < size; ++i ) {
    int j = i - 1;
    int k = i + 1;

    if ( i == 0 ) {
      averageNNdistance += fabs( *( evts[i] + index ) - *( evts[k] + index ) );
      continue;
    }
    if ( k == size ) {
      averageNNdistance += fabs( *( evts[i] + index ) - *( evts[j] + index ) );
      continue;
    } else {
      double low  = fabs( *( evts[i] + index ) - *( evts[k] + index ) );
      double high = fabs( *( evts[i] + index ) - *( evts[j] + index ) );
      averageNNdistance += low > high ? high : low;
    }
  }
  double avg = averageNNdistance / double( size );

  double sigma = 0;
  for ( int i = 0; i < size; ++i ) {
    int j = i - 1;
    int k = i + 1;
    if ( i == 0 ) {
      sigma += pow( fabs( *( evts[i] + index ) - *( evts[k] + index ) ) - avg, 2 );
      continue;
    }
    if ( k == size ) {
      sigma += pow( fabs( *( evts[i] + index ) - *( evts[j] + index ) ) - avg, 2 );
      continue;
    } else {
      double low  = fabs( *( evts[i] + index ) - *( evts[k] + index ) );
      double high = fabs( *( evts[i] + index ) - *( evts[j] + index ) );
      sigma += fabs( low ) > fabs( high ) ? pow( high - avg, 2 ) : pow( low - avg, 2 );
    }
  }
  return sigma / double( evts.size() );
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

  auto sorter = [&index]( double* a, double* b ) { return *( a + index ) < *( b + index ); };
  DEBUG( "Sorting" );
  std::sort( evts.begin(), evts.end(), sorter );

  unsigned int midpoint = evts.size() / 2;
  DEBUG( "Done sorting, midpoint = " << midpoint );
  double midposition =
    evts.size() % 2 == 0
    ? ( *( evts[midpoint - 1] + index ) + *( evts[midpoint] + index ) ) / 2
    : ( *( evts[midpoint + 1] + index ) + *( evts[midpoint - 1] + index ) + *( evts[midpoint] + index ) ) / 3;
  DEBUG( "Splitting on midpoint = " << midpoint );
  midpoint += midposition > *( evts[midpoint] + index );

  std::vector<double*> leftEvents( evts.begin(), evts.begin() + midpoint );
  std::vector<double*> rightEvents( evts.begin() + midpoint, evts.end() );
  DEBUG( "Split list into: " << leftEvents.size() << " " << rightEvents.size() );
  indexQueue.pop();
  if ( indexQueue.empty() ) refreshQueue( evts, indexQueue, depth + 1 );
  auto left  = makeNodes( leftEvents, indexQueue, depth + 1 );
  auto right = makeNodes( rightEvents, indexQueue, depth + 1 );
  auto node  = std::make_shared<Decision>( index, midposition, right, left );
  return node;
}

std::shared_ptr<BinDT::INode> BinDT::makeNodes( std::vector<double*> source, std::vector<double*> target )
{
  std::queue<unsigned int> iq;
  refreshQueue( source, iq, 0 );
  auto node = makeNodes( source, target, iq, 0 );
  m_top     = node;
  return node;
}

double BinDT::getBestPost( std::vector<double*> source, std::vector<double*> target, int index, bool verbose )
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
  if ( verbose ) {
    INFO( "Checking: " << source.size() << " " << target.size() );
  }
  for ( unsigned int i = 0; i < target.size() - m_minEvents; i++ ) {
    targetLeftWeight += target[i][m_dim];
    targetRightWeight -= target[i][m_dim];
    if ( i < m_minEvents ) continue;
    double pos = ( target[i][index] + target[i + 1][index] ) / 2.;

    for ( ; sourceIt != source.end(); ++sourceIt ) {
      double positionOfThis = *( *sourceIt + index );
      if ( verbose )
        INFO( "Position = " << pos << " current = " << positionOfThis << " population = " << sourceLeftPopulation );
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
  std::sort( source.begin(), source.end(), sorter );
  std::sort( target.begin(), target.end(), sorter );
  double optPos = getBestPost( source, target, index );

  std::vector<double*> sourceLeft;  //( source.begin(), optSourceIt );
  std::vector<double*> targetLeft;  //( target.begin(), optTargetIt );
  std::vector<double*> sourceRight; //( optSourceIt, source.end() );
  std::vector<double*> targetRight; //( optTargetIt, target.end() );
  for ( auto& event : source ) {
    ( event[index] > optPos ? sourceRight : sourceLeft ).push_back( event );
  }
  for ( auto& event : target ) {
    ( event[index] > optPos ? targetRight : targetLeft ).push_back( event );
  }
  if ( sourceLeft.size() == 0 || sourceRight.size() == 0 ) {
    DEBUG( "No data to divide into block" );
    DEBUG( "Dividing:  source=(" << sourceLeft.size() << ", " << sourceRight.size() << "); target=("
        << targetLeft.size() << ", " << targetRight.size() << "); index=" << index
        << " depth = " << depth );
    auto node = std::make_shared<EndNode>( m_endNodes.size() );
    m_endNodes.push_back( node );
    return node;
    /// getBestPost( source, target, index, true);
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
    for ( unsigned int i = 0; i < m_dim; ++i ) indices.emplace_back( i, nnUniformity( evts, i ) );
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
