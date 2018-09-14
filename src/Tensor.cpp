#include "AmpGen/Tensor.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>

#include "AmpGen/Utilities.h"

using namespace AmpGen;

std::ostream& AmpGen::operator <<( std::ostream& out, const Tensor::Index& index ){
  return out << index.m_ptr << " " << index.m_isUpper; 
}  
Tensor::Tensor()
{
  DEBUG("Begin constructor <default>");
  setupCoordinates();
  DEBUG("End constructor <default>");
}

Tensor::Tensor( const std::vector<size_t>& _dim ) 
  : m_dim( _dim ), m_elements( nElements(), Constant( 0. ) )
{
  DEBUG("Begin constructor <empty>");
  setupCoordinates();
  DEBUG("End constructor <empty>");
}

Tensor::Tensor( const std::vector<Expression>& elements )
  : m_dim( std::vector<size_t>( {elements.size()} ) )
{
  DEBUG("Begin constructor <Expressions>");
  setupCoordinates();
  for( auto& element : elements ) append( element );
  DEBUG("End constructor <Expressions>");
}

Expression Tensor::get( const size_t& co )
{
  if ( co >= m_elements.size() )
    ERROR( "Element (" + std::to_string( co ) + " ) out of range (0"
        << ", " << m_elements.size() << ")" );
  return ( m_elements[m_symmetrisedCoordinates[co]] );
}
Expression Tensor::get( const size_t& co ) const
{
  if ( co >= m_elements.size() )
    ERROR( "Element (" + std::to_string( co ) + " ) out of range (0"
        << ", " << m_elements.size() << ")" );
  return ( m_elements[m_symmetrisedCoordinates[co]] );
}

std::string Tensor::to_string(const ASTResolver* resolver) const {
  std::string value = "{";
  for(size_t i = 0 ; i < size(); ++i )
  {
    value += Tensor::operator[](i).to_string(resolver) + (i == size() -1  ? "}" : ", " ) ; 
  }
  return value;
}

size_t Tensor::rank() const { return m_dim.size(); }

int Tensor::metricSgn( const std::vector<size_t>& coordinates ) const
{
  int sgn = 1;
  for ( auto& coord : coordinates ) sgn *= ( coord == 3 ) ? 1 : -1;
  return sgn;
}

int Tensor::metricSgn( const size_t& index ) const { return metricSgn( coords( index ) ); }

void Tensor::append( const Expression& expression )           { m_elements.emplace_back( expression ); }
void Tensor::append( const Parameter& parameter )             { m_elements.emplace_back( parameter ); }
void Tensor::append( const double& value )                    { m_elements.emplace_back( Constant( value )); }
void Tensor::append( const std::complex<double>& value )      { m_elements.emplace_back( Constant( value )); }
void Tensor::append( const std::string& name, bool resolved ) { m_elements.emplace_back( Parameter( name, 0, resolved )); }

void Tensor::setupCoordinates()
{
  int st=1;
  for( auto& d : m_dim ) st *= d; 
  DEBUG( "Setting up coordinates: " << dimString() << " coordinate mapping = " << st );
  m_symmetrisedCoordinates.resize( st  );
  std::iota( m_symmetrisedCoordinates.begin(), m_symmetrisedCoordinates.end(), 0 );
  m_uniqueElements.resize( st ) ; 
  std::iota( m_uniqueElements.begin(), m_uniqueElements.end(), 0 );
}

Expression Tensor::get( const std::vector<size_t>& _co ) const { return ( m_elements[index( _co )] ); }

size_t Tensor::size() const { return m_elements.size(); }

size_t Tensor::index( const std::vector<size_t>& _co ) const
{
  auto id = Tensor::coordinates_to_index( _co, m_dim );
  if ( id > nElements() ) ERROR( "Element (" << id << ") out of range" );
  return id;
}

size_t Tensor::symmetrisedIndex( const std::vector<size_t>& _co ) const 
{
  auto id = Tensor::coordinates_to_index( _co, m_dim );
  return m_symmetrisedCoordinates[ id ] ;
}

const std::vector<size_t> Tensor::coords( const size_t& index ) const
{
  return Tensor::index_to_coordinates( index, m_dim );
}

std::vector<size_t> Tensor::index_to_coordinates( const size_t& index, const std::vector<size_t>& dim )
{
  std::vector<size_t> returnValue;
  size_t index_temp = index;
  for (size_t j = 1; j < dim.size() + 1; ++j ) {
    size_t dproduct = 1;
    for (size_t i = 0; i < dim.size() - j; ++i ) {
      dproduct *= dim[i];
    }
    size_t val = ( index_temp - ( index_temp % dproduct ) ) / dproduct;
    index_temp -= dproduct * val;
    returnValue.push_back( val );
  }
  return returnValue;
}

std::string Tensor::coordinates_to_string( const std::vector<size_t>& coordinates )
{
  std::string rt = "[";
  for(size_t i = 0; i < coordinates.size(); ++i)
    rt += std::to_string( coordinates[i] )+ ( i==coordinates.size()-1?"]":",");
  return rt;
}


size_t Tensor::coordinates_to_index( const std::vector<size_t>& _co, const std::vector<size_t>& dim )
{
  size_t _index   = 0;
  size_t dproduct = 1;
  for ( size_t i = 0; i != dim.size() ; ++i ) {
    _index += _co[dim.size() -1 -i ] * dproduct;
    dproduct *= dim[dim.size()-i-1];
  }
  return _index;
}

size_t Tensor::nElements() const
{
  size_t dim = 1;
  for ( auto& d : m_dim ) {
    dim *= ( d != 0 ) ? d : 1;
  }
  return dim;
}
const std::string Tensor::dimString() const 
{
  std::string str="[";
  for( auto it : m_dim ) str += std::to_string(it) + ", ";
    str = str.substr(0, str.size() -2 );
  return str + "]";
}

size_t Tensor::nDim() const { return m_dim.size(); }

bool Tensor::rankMatches( const Tensor& other )
{
  bool success = true;
  if ( m_dim.size() != other.m_dim.size() ) return false;
  for ( size_t i = 0; i < m_dim.size(); ++i ) success &= m_dim[i] == other.m_dim[i];
  return success;
}

Tensor AmpGen::operator+( const Tensor& t1, const Tensor& t2 )
{
  Tensor result( t1.dims() );

  if ( t1.rank() != t2.rank() ) {
    ERROR( "Addition between tensors of different rank makes no sense!" );
    t1.print();
    t2.print();
  }
  if ( t1.nElements() != t2.nElements() ) {
    ERROR( "Addition between tensors with different number of elements " << t1.nElements() << " " << t2.nElements() );
  }

  for ( size_t i = 0; i < t1.nElements(); ++i ) result[i] = t1[i] + t2[i];

  return result;
}

Tensor AmpGen::operator-( const Tensor& t1, const Tensor& t2 )
{
  Tensor result( t1.dims() );
  if ( t1.rank() != t2.rank() ) {
    ERROR( "Subtraction between tensors of different rank makes no sense!" );
    t1.print();
    t2.print();
  }
  if ( t1.nElements() != t2.nElements() ) {
    t1.print();
    t2.print();
    ERROR( "Subtraction between tensors with different number of elements " << t1.nElements() << " "
        << t2.nElements() );
  }
  for ( size_t i = 0; i < t1.nElements(); ++i ) result[i] = t1[i] - t2[i];
  if ( result.nElements() != t1.nElements() || result.nElements() != t2.nElements() ) {
    ERROR( "What are you saying?" );
  }

  return result;
}

void Tensor::st(const bool simplify) {
  for(size_t i = 0 ; i < size(); ++i ){
    m_elements[i] = make_cse( m_elements[i], simplify ); 
  }
}

Tensor Tensor::conjugate() const 
{
  Tensor copy( dims() );
  for(size_t i = 0 ; i < size(); ++i ) copy[i] = fcn::conj( get(i) );
  return copy; 
}

Tensor AmpGen::operator/( const Tensor& t1, const Expression& t2 )
{
  Tensor result( t1.dims() );
  for (size_t i = 0; i < t1.nElements(); ++i ) result[i] = t1[i] / t2;

  return result;
}

Tensor AmpGen::operator*( const Expression& other, const Tensor& t1 )
{
  Tensor result( t1.dims() );
  for (size_t i = 0; i < t1.nElements(); ++i ) result[i] = t1[i] * other;
  return result;
}

Tensor AmpGen::operator*( const Tensor& t1, const Expression& other ) { return other * t1; }

Tensor AmpGen::operator/( const Tensor& t1, const double& t2 ) { return t1 / Constant( t2 ); }
Tensor AmpGen::operator*( const double& other, const Tensor& t1 ) { return Constant( other ) * t1; }
Tensor AmpGen::operator*( const Tensor& t1, const double& other ) { return t1 * Constant( other ); }

void AmpGen::Tensor::operator+=( const Tensor& rhs ){ *this = *this + rhs; }
void AmpGen::Tensor::operator-=( const Tensor& rhs ){ *this = *this - rhs; }
Expression AmpGen::dot( const Tensor& A, const Tensor& B )
{

  if ( A.nElements() == 1 && B.nElements() == 1 ) return A[0] * B[0];

  if ( A.rank() != B.rank() || A.nElements() != B.nElements() ) {
    ERROR( "Rank(A) = " << A.rank() << " != Rank(B) = " << B.rank() );
    ERROR( "Elem(A) = " << A.nElements() << " ; Elem(B) = " << B.nElements() );
    ERROR( "Cannot fully contract unmatched Tensors - you nitwit" );
    return Constant( 0 );
  }
  Expression result;
  for ( int i = A.nElements() -1 ; i != -1 ; --i ) {
    int sgn = A.metricSgn( i );
    if ( sgn == 1 )
      result = result + A.get( i ) * B.get( i );
    else if ( sgn == -1 )
      result = result - A.get( i ) * B.get( i );
    else
      result + sgn* A.get( i ) * B.get( i );
  }
  return result;
}

Tensor Tensor::Invert() const
{
  int actualsize = m_dim[0];
  if ( rank() != 2 ) {
    ERROR( "Inversion only implemented for rank 2 objects" );
  };
  Tensor data = *this;
  for ( int i = 1; i < actualsize; i++ ) data[i] = data[i] / data[0]; // normalize row 0
  // data.print();

  for ( int i = 1; i < actualsize; i++ ) {
    for ( int j = i; j < actualsize; j++ ) { // do a column of L
      Expression sum = 0;
      for ( int k = 0; k < i; k++ ) sum = sum + data[j * actualsize + k] * data[k * actualsize + i];
      //   INFO("Subtracting " << sum.to_string() << "   " << data[j*actualsize+i].to_string() );
      data[j * actualsize + i] = data[j * actualsize + i] - SubTree( sum );
    }
    if ( i == actualsize - 1 ) continue;
    for ( int j = i + 1; j < actualsize; j++ ) { // do a row of U
      Expression sum = 0.0;
      for ( int k = 0; k < i; k++ ) sum = sum + data[i * actualsize + k] * data[k * actualsize + j];
      data[i * actualsize + j] = ( data[i * actualsize + j] - SubTree( sum ) ) / data[i * actualsize + i];
    }
  }
  // INFO("LU part : " ) ;
  // data.print();
  for ( int i = 0; i < actualsize; i++ ) // invert L
    for ( int j = i; j < actualsize; j++ ) {
      Expression x = 1.0;
      if ( i != j ) {
        x = 0.0;
        for ( int k = i; k < j; k++ ) x = x - data[j * actualsize + k] * data[k * actualsize + i];
      }
      data[j * actualsize + i] = x / data[j * actualsize + j];
    }
  for ( int i = 0; i < actualsize; i++ ) // invert U
    for ( int j = i; j < actualsize; j++ ) {
      if ( i == j ) continue;
      Expression sum = 0.0;
      for ( int k = i; k < j; k++ )
        sum = sum + data[k * actualsize + j] * ( ( i == k ) ? 1.0 : data[i * actualsize + k] );
      data[i * actualsize + j] = -sum;
    }
  for ( int i = 0; i < actualsize; i++ ) // final inversion
    for ( int j = 0; j < actualsize; j++ ) {
      Expression sum = 0.0;
      for ( int k = ( ( i > j ) ? i : j ); k < actualsize; k++ )
        sum = sum + ( ( j == k ) ? 1.0 : data[j * actualsize + k] ) * data[k * actualsize + i];
      data[j * actualsize + i] = sum;
    }
  return data;
}

/// generic contraction operator ///

void Tensor::print() const
{
  std::string dim_string = "";
  for ( unsigned int i = 0; i < m_dim.size(); ++i )
    dim_string += std::to_string( m_dim[i] ) + ( i == m_dim.size() - 1 ? "" : "x" );
  INFO( "Rank=" << m_dim.size() << ", Dim=(" << dim_string << ")" );
  if ( m_dim.size() == 0 )
    std::cout << m_elements[0].to_string() << std::endl;
  else if ( m_dim.size() == 1 ) {
    for ( unsigned int x = 0; x < m_dim[0]; ++x ) std::cout << m_elements[x].to_string() << std::endl;
  }

  else if ( m_dim.size() == 2 ) {
    for ( unsigned int y = 0; y < m_dim[1]; ++y ) {
      for ( unsigned int x = 0; x < m_dim[0]; ++x ) {
        std::cout << m_elements[index( {x, y} )].to_string() << "     ";
      }
      std::cout << std::endl;
    }
  } else {
    for ( unsigned int i = 0; i < m_elements.size(); ++i ) {
      INFO( "M(" + vectorToString( coords( i ) ) + ") = " << m_elements[i].to_string() );
    }
  }
}

TensorProxy Tensor::operator()( const std::vector<Tensor::Index>& indices ) const
{
  return TensorProxy( *this, indices );
}
TensorProxy Tensor::operator()( const Tensor::Index& a ) const { return TensorProxy( *this, {a} ); }
TensorProxy Tensor::operator()( const Tensor::Index& a, const Tensor::Index& b ) const
{
  return TensorProxy( *this, {a, b} );
}
TensorProxy Tensor::operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c ) const
{
  return TensorProxy( *this, {a, b, c} );
}

TensorProxy Tensor::operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c,
    const Tensor::Index& d ) const
{
  return TensorProxy( *this, {a, b, c, d} );
}

TensorProxy AmpGen::operator*( const TensorProxy& t1, const TensorProxy& t2 )
{
  struct contractor {
    unsigned int i;
    unsigned int j;
    int sgn;
    contractor( const size_t& i, const size_t&j, const int& sgn) : i(i),j(j),sgn(sgn) {}
  };
  std::vector<contractor> contractions;

  std::vector<Tensor::Index> unsummedIndices;

  const std::vector<Tensor::Index>& t1_index = t1.indices();
  const std::vector<Tensor::Index>& t2_index = t2.indices();
  const size_t t1_size                       = t1.indices().size();
  const size_t t2_size                       = t2.indices().size();

  const Tensor& t1_tensor = t1.tensor();
  const Tensor& t2_tensor = t2.tensor();

  std::vector<size_t> finalTensorRank;
  std::vector<size_t> contractionMatrix;

  for( size_t i = 0; i < t1_size; ++i )
  {
    if ( std::find( t2_index.begin(), t2_index.end(), t1.indices()[i] ) == t2_index.end() ) {
      unsummedIndices.push_back( t1.indices()[i] );
      finalTensorRank.push_back( t1_tensor.dims()[i] );
    }
  }
  for( size_t i = 0; i < t2_size; ++i ) {
    auto it = std::find( t1_index.begin(), t1_index.end(), t2_index[i] );
    if ( it == t1_index.end() ) {
      unsummedIndices.push_back( t2.indices()[i] );
      finalTensorRank.push_back( t2_tensor.dims()[i] );
    } else {
      contractions.emplace_back( 
      std::distance( t1_index.begin(), it ), i, 
      t2_index[i].isUpper() != it->isUpper() ? -1 : 1 );
    }
  }
  unsigned int sumElements = 1;
  for ( auto& c : contractions ) {
    unsigned int dim = t1_tensor.dims()[c.i];
    contractionMatrix.push_back( dim );
    sumElements *= dim;
  }
  Tensor value( finalTensorRank );

  unsigned nElem = value.nElements();
  DEBUG( "Got " << t1_tensor.dims().size() << " x " << t2_tensor.dims().size() << " with " << contractions.size() << " contractions " );

  for ( unsigned int elem = 0; elem < nElem; ++elem ) {
    auto coords = Tensor::index_to_coordinates( elem, finalTensorRank ); // value.coords(elem); /// coordinates of this
    // part of the expression ///
    std::vector<size_t> t1_coords( t1_size, 0 );
    std::vector<size_t> t2_coords( t2_size, 0 );
    unsigned int i = 0;
    unsigned int j = 0;
    do {
      //DEBUG( "Checking : " << i );
      if ( !isIn( contractions, i, [](
              // INFO( "checking : " << i );
              const contractor& a, const unsigned int& b ) { return a.i == b; } ) ) {
      //  DEBUG( "Setting coordinate " << i << " = " << j << " coords.size() = " << coords.size() );
        t1_coords[i] = coords[j];
        j++;
      }
    } while ( ++i < t1_size );
    i = 0;
    j = t1_size - contractions.size();
    do {
     // DEBUG( "Checking : " << i );
      if ( !isIn( contractions, i, []( const contractor& a, const unsigned int& b ) { return a.j == b; } ) ) {
        t2_coords[i] = coords[j];
      //  DEBUG( "Setting i = " << i << ", (" << j << ")" << coords[j] );
        j++;
      }
    } while ( ++i < t2_size );
    /// calculate the number of components to sum over.
    Expression elementExpression = 0;
    for ( unsigned int i = 0; i < sumElements; ++i ) {
      auto contractedCoordinates = Tensor::index_to_coordinates( i, contractionMatrix );
      int sign                   = 1;
      for ( unsigned int i = 0; i < contractions.size(); ++i ) {
        t1_coords[contractions[i].i] = contractedCoordinates[i];
        t2_coords[contractions[i].j] = contractedCoordinates[i];
        sign *= ( contractedCoordinates[i] ) == 3 ? 1 : contractions[i].sgn;
      }
      elementExpression = elementExpression + sign * t1_tensor[t1_coords] * t2_tensor[t2_coords];
    }
    value[elem] = elementExpression;
  }
//  DEBUG( "Returning object of rank = " << unsummedIndices.size() );
  return TensorProxy( value, unsummedIndices );
}

TensorProxy AmpGen::operator/( const TensorProxy& t1, const Expression& t2 )
{
  return TensorProxy( t1.tensor() / t2, t1.indices() );
}
TensorProxy AmpGen::operator*( const Expression& t1, const TensorProxy& t2 )
{
  return TensorProxy( t1 * t2.tensor(), t2.indices() );
}
TensorProxy AmpGen::operator*( const TensorProxy& t1, const Expression& t2 )
{
  return TensorProxy( t1.tensor() * t2, t1.indices() );
}

TensorProxy AmpGen::operator*( const TensorProxy& t1, const double& t2 )
{
  return TensorProxy( t1.tensor() * Constant( t2 ), t1.indices() );
}
TensorProxy AmpGen::operator*( const double& t2, const TensorProxy& t1 )
{
  return TensorProxy( t1.tensor() * Constant( t2 ), t1.indices() );
}

TensorProxy AmpGen::operator/( const TensorProxy& t1, const double& t2 )
{
  return TensorProxy( t1.tensor() / Constant( t2 ), t1.indices() );
}

TensorProxy AmpGen::operator-( const TensorProxy& t1, const TensorProxy& t2 ) { return t1 + ( -1 ) * t2; }

std::string coordinateString( const std::vector<size_t>& indices ){
  std::string rt = "";
  for( auto& index : indices ) rt += std::to_string(index) + ", ";
  rt.erase( rt.end()-2, rt.end() );
  return rt ;
}

TensorProxy AmpGen::operator+( const TensorProxy& t1, const TensorProxy& t2 )
{
  auto indices_t1 = t1.indices();
  auto indices_t2 = t2.indices();
  auto tensor_t1  = t1.tensor();
  auto tensor_t2  = t2.tensor();
  std::vector<unsigned int> t1_to_t2_mapping;
  for ( unsigned int j = 0; j < indices_t2.size(); ++j ) {
    for ( unsigned int i = 0; i < indices_t1.size(); ++i ) {
      if ( indices_t1[i] == indices_t2[j] ){
        t1_to_t2_mapping.push_back( i );
        break;
      }
    }
  }
  if ( t1_to_t2_mapping.size() != indices_t1.size() ) {
    ERROR( "Mapping illformed!" << indices_t1.size() << " " << indices_t2.size()  );
    return TensorProxy( Tensor(), {} );
  }
  Tensor result( tensor_t1 );

  for ( auto& i : tensor_t1.uniqueElements() ) {
    auto t1_coords = tensor_t1.coords( i );
    auto t2_coords = std::vector<size_t>( t1_coords.size(), 0 );
    for( unsigned int j=0;j<t1_coords.size(); ++j ) t2_coords[ j ] = t1_coords[t1_to_t2_mapping[j]];
    result[i] = result[i] + tensor_t2[t2_coords];
  }
  return TensorProxy( result, indices_t1 );
}

void Tensor::imposeSymmetry( size_t indexA, size_t indexB) 
{
  if( indexA > indexB ) std::swap( indexA, indexB );
  std::map< size_t, size_t > counter; 
  for( size_t i = 0 ; i < m_symmetrisedCoordinates.size(); ++i ){
    auto coordinates = Tensor::index_to_coordinates( i, m_dim ); /// raw coordinates of this ///
    if( coordinates[indexB] > coordinates[indexA] ) 
      std::swap( coordinates[indexA], coordinates[indexB] );
    auto symmetrisedIndex = Tensor::coordinates_to_index( coordinates, m_dim );
    counter[m_symmetrisedCoordinates[symmetrisedIndex]]++;
    m_symmetrisedCoordinates[i] = m_symmetrisedCoordinates[symmetrisedIndex];
  }
  m_uniqueElements.clear();
  for( auto& f : counter ) m_uniqueElements.push_back( f.first );
  DEBUG("Imposing symmetries on " << indexA << " <--> " << indexB << " reduces size to: " << counter.size() );
}

void Tensor::imposeSymmetry( std::vector<size_t> indices )
{
  std::sort( indices.begin(), indices.end() );
  for( size_t i=0;i<indices.size()-1;++i){
    for( size_t j=i+1;j<indices.size();++j){
      imposeSymmetry( indices[i], indices[j] );
    }
  }
}

TensorProxy::TensorProxy( const Tensor& tensor, const std::vector<Tensor::Index>& indices ) : m_tensor( tensor )
{
  if ( m_tensor.rank() != indices.size() ) {
    ERROR( "Wrong number of indices set (#Indices=" << indices.size() << ", rank=" << m_tensor.nDim() << ")");
    return;  
  } 
  m_indices = indices; 
  std::vector< std::pair<size_t,size_t> > contractions; 
  for( size_t i = 0 ; i < indices.size(); ++i ){
    size_t j = i + 1;
    for( ; j < indices.size(); ++j ){
      if( indices[i] == indices[j] ){
        WARNING("Trace-like contractions on Tensors not implemented");
        contractions.emplace_back( i, j ); 
        break ; 
      } 
    }
   // if( j == indices.size() ) m_indices.push_back( indices[i] );
  }
  
}

TensorProxy TensorProxy::reorder( const std::vector<Tensor::Index>& indices )
{
  std::vector<size_t> mapping( m_indices.size() ,0 );
  for(size_t j= 0; j < indices.size(); ++j)
  {
    for(size_t i = 0; i < m_indices.size(); ++i)
    {
      if( m_indices[j] == indices[i] ) mapping[j] = i;
    }
  }
  Tensor reordered( m_tensor.dims() ); 
  for(size_t i = 0 ; i < m_tensor.size(); ++i)
  { 
    auto coordinates       = Tensor::index_to_coordinates(i,m_tensor.dims() );
    std::vector<size_t> new_coordinates( m_indices.size() );
    for(size_t z = 0 ; z < m_indices.size(); ++z )
      new_coordinates[z] = coordinates[mapping[z]];
    reordered[ i ] = m_tensor[ new_coordinates ];
  }
  return TensorProxy( reordered, indices );
}  

Tensor AmpGen::Identity( const size_t& rank )
{
  Tensor id( std::vector<size_t>( {rank, rank} ) );
  for (size_t i = 0; i < rank; ++i) id(i, i) = 1;
  return id;
}

const Tensor AmpGen::LeviCivita( const size_t& rank )
{
  std::vector<size_t> dim( rank, rank );
  std::vector<size_t> indices( rank );

  std::iota( indices.begin(), indices.end(), 0 );
  auto permutation_sign = []( const std::vector<size_t>& permutation){
    int product = 1;
    for ( size_t i = 0; i < permutation.size() - 1; ++i ) {
      for ( size_t j = i + 1; j < permutation.size(); ++j ) {
        product *= ( (int)permutation[i] - (int)permutation[j] );
      }
    }
    return product;
  };
  int p0 = permutation_sign( indices );

  Tensor result( dim ); /// create tensor of rank N ///
  do {
    size_t index = result.index( indices );
    result[index]      = permutation_sign( indices ) / p0;
  } while ( std::next_permutation( indices.begin(), indices.end() ) );
  return result;
}

