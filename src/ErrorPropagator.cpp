#include "AmpGen/ErrorPropagator.h"

#include <ostream>

#include "AmpGen/MinuitParameterSet.h"
#include "TDecompChol.h"
#include "TRandom3.h"

using namespace AmpGen;

GaussErrorPropagator::GaussErrorPropagator( const TMatrixD& reducedCovariance, const std::vector<MinuitParameter*>& params, TRandom3* rnd )
  : m_parameters( params ), m_rand( rnd ), m_decomposedCholesky( params.size(), params.size() )
{
  for ( size_t x = 0; x < params.size(); ++x ) {
    auto p = params[x];
    INFO( p->name() << "  " << p->mean() << " +/- " << sqrt( reducedCovariance( x, x ) ) );
    m_startingValues.push_back( p->mean() );
  }
  TDecompChol decomposed( reducedCovariance );
  decomposed.Decompose();
  m_decomposedCholesky = decomposed.GetU();
  /// transpose the cholesky matrix
  for ( int i = 0; i < m_decomposedCholesky.GetNrows(); ++i ) {
    for ( int j = i + 1; j < m_decomposedCholesky.GetNrows(); ++j ){
      std::swap( m_decomposedCholesky(i, j), m_decomposedCholesky(j, i)  );
    }
  }
}

void GaussErrorPropagator::perturb()
{
  const unsigned int N = m_decomposedCholesky.GetNrows();
  TVectorD e( N );
  for ( unsigned int i = 0; i < N; ++i ) e[i] = m_rand->Gaus( 0, 1 );
  TVectorD p = m_decomposedCholesky * e; 
  for ( int j = 0; j < p.GetNrows(); ++j ) {
    auto f = m_parameters[j];
    f->setCurrentFitVal( m_startingValues[j] + p[j] );
  }
}

void GaussErrorPropagator::reset()
{
  for ( unsigned int j = 0; j < m_parameters.size(); ++j ) m_parameters[j]->setCurrentFitVal( m_startingValues[j] );
}

LinearErrorPropagator::LinearErrorPropagator( const TMatrixD& reducedCovarianceMatrix,
    const std::vector<MinuitParameter*>& params )
  : m_cov( reducedCovarianceMatrix ), m_parameters( params )
{
}

LinearErrorPropagator::LinearErrorPropagator( const std::vector<MinuitParameter*>& params )
{
  for( auto& param : params ){
    if( !param->isFree() || param->err() == 0 ) continue;
    m_parameters.push_back( param );
  }
  m_cov.ResizeTo( m_parameters.size(), m_parameters.size() );
  for( size_t i = 0 ; i < m_parameters.size(); ++i ) 
    m_cov(i,i) = m_parameters[i]->err() * m_parameters[i]->err();
}

LinearErrorPropagator::LinearErrorPropagator( const MinuitParameterSet& mps )
{
  for(auto& param : mps){
    if( ! param->isFree() || param->err() == 0 ) continue; 
    m_parameters.push_back(param);
  }
  m_cov.ResizeTo( m_parameters.size(), m_parameters.size() );
  for( size_t i = 0 ; i < m_parameters.size(); ++i ) 
    m_cov(i,i) = m_parameters[i]->err() * m_parameters[i]->err();
}

void LinearErrorPropagator::add( const LinearErrorPropagator& p2 )
{
  size_t superSet = size();
  size_t oldSize  = size();
  auto p1_pMap    = posMap();
  std::vector<size_t> props( p2.size() );
  for ( size_t x = 0; x != p2.size(); ++x ) {
    auto it = p1_pMap.find( p2.params()[x]->name() );
    if ( it == p1_pMap.end() ) {
      props[x] = superSet++;
      m_parameters.push_back( p2.params()[x] );
    } else
      props[x] = it->second;
  }
  TMatrixD old_cov = m_cov;
  if ( superSet != oldSize ) m_cov.ResizeTo( superSet, superSet );
  for ( size_t x = 0; x != oldSize; ++x ) {
    for ( size_t y = 0; y != oldSize; ++y ) m_cov( x, y ) = old_cov( x, y );
  }
  auto p2_cov = p2.cov();
  for ( size_t x = 0; x < p2.size(); ++x ) {
    for ( size_t y = 0; y < p2.size(); ++y ) {
      auto xp = props[x];
      auto yp = props[y];
      m_cov( xp, yp ) = m_cov( xp, yp ) + p2_cov( x, y );
    }
  }
}

double LinearErrorPropagator::getError( const std::function<double(void)>& fcn ) const
{
  unsigned int N = m_cov.GetNrows();
  TVectorD errorVec( N );
  for ( unsigned int i = 0; i < N; ++i ) {
    DEBUG( "Perturbing parameter: [" << m_parameters[i]->name() << "] " << m_parameters[i] << " by "
        << sqrt( m_cov( i, i ) ) << " " << m_parameters[i] );
    errorVec(i) = derivative(fcn,i);
    fcn(); 
  }
  return sqrt( errorVec * ( m_cov * errorVec ) );
}

std::vector<double> LinearErrorPropagator::getVectorError( const std::function<std::vector<double>(void)>& fcn, size_t RANK ) const
{
  unsigned int N = m_cov.GetNrows();
  std::vector<TVectorD> errorVec( RANK, TVectorD( N ) );
  for ( unsigned int i = 0; i < N; ++i ) {
    double startingValue = m_parameters[i]->mean();
    double error         = sqrt( m_cov( i, i ) );
    double min           = m_parameters[i]->mean() - error;
    double max           = m_parameters[i]->mean() + error;
    DEBUG( "Perturbing parameter: " << m_parameters[i]->name() << " -> [" << min << ", " << max << "]" );

    m_parameters[i]->setCurrentFitVal( max );
    auto plus_variation = fcn();
    m_parameters[i]->setCurrentFitVal( min );
    auto minus_variation = fcn();
    m_parameters[i]->setCurrentFitVal( startingValue );
    for ( size_t j = 0; j < RANK; ++j ) {
      errorVec[j]( i ) = ( plus_variation[j] - minus_variation[j] ) / ( 2 * error );
    }
  }
  fcn();
  std::vector<double> rt( RANK, 0 );
  for ( unsigned int j = 0; j < RANK; ++j ) rt[j] = sqrt( errorVec[j] * ( m_cov * errorVec[j] ) );
  return rt;
}


void LinearErrorPropagator::reset()
{
  for ( int x = 0; x < m_cov.GetNrows(); ++x ) {
    for ( int y = 0; y < m_cov.GetNrows(); ++y ) {
      m_cov( x, y ) = 0;
    }
  }
  m_parameters.clear();
}

TMatrixD LinearErrorPropagator::covarianceMatrix(const std::vector<std::function<double(void)>>& functions ){
  size_t M = functions.size();
  size_t N = size();
  TMatrixD A(M,N);
  for( size_t k = 0 ; k < M; ++k )
    for( size_t i = 0 ; i < N; ++i ) 
      A(k,i) = derivative( functions[k], i );

  TMatrixD vci( M,M);

  for( size_t i = 0; i < M ; ++i ){
    for( size_t j = 0; j < M ; ++j ){
      for(size_t k = 0 ; k < N ; ++k ){
        for(size_t l = 0 ; l < N ; ++l ){
          vci(i,j) += A(i,k) * m_cov(k,l) * A(j,l);
        }
      }
    }
  }
  return vci;
}

std::pair<double, double> LinearErrorPropagator::combinationCovWeighted( 
    const std::vector<std::function<double(void)>>& functions)
{
  auto cov = covarianceMatrix( functions );
  TMatrixD covInverse = cov; 
  covInverse.Invert();
  double wt     = 0;
  double pred   = 0;
  double sigma2 = 0;
  for( size_t i = 0 ; i < functions.size(); ++i ){
    for( size_t j = 0 ; j < functions.size(); ++j ){
      wt += covInverse(i,j);
    }
  }
  std::vector<double> weights( functions.size() );
  for( size_t i = 0 ; i < functions.size(); ++i ){
    for( size_t j = 0 ; j < functions.size(); ++j ) weights[i] += covInverse(i,j);
    pred += functions[i]() * weights[i] / wt; 
  }
  
  for( size_t i = 0 ; i < functions.size(); ++i ){
    for( size_t j = 0 ; j < functions.size(); ++j ){
      sigma2 += weights[i] * weights[j] * cov(i,j) / (wt*wt);
    }
  }
  return std::make_pair( pred, sqrt(sigma2) );
}

std::vector<double> LinearErrorPropagator::combinationWeights( const std::vector< std::function<double(void)>>& estimators )
{
  auto cov = covarianceMatrix( estimators );
  TMatrixD covInverse = cov; 
  covInverse.Invert();
  double wt     = 0;
  for( size_t i = 0 ; i < estimators.size(); ++i ){
    for( size_t j = 0 ; j < estimators.size(); ++j ){
      wt += covInverse(i,j);
    }
  }
  std::vector<double> weights( estimators.size() );
  for( size_t i = 0 ; i < estimators.size(); ++i ){
    for( size_t j = 0 ; j < estimators.size(); ++j ) weights[i] += covInverse(i,j)/wt;
  }
  return weights; 
}

size_t LinearErrorPropagator::size() const { return m_parameters.size(); }
const TMatrixD& LinearErrorPropagator::cov() const { return m_cov; }
const std::map<std::string, size_t> LinearErrorPropagator::posMap() const
{
  std::map<std::string, size_t> pMap;
  for ( size_t pos = 0; pos != m_parameters.size(); ++pos ) {
    pMap[m_parameters[pos]->name()] = pos;
  }
  return pMap;
}

TMatrixD LinearErrorPropagator::correlationMatrix(const std::vector<std::function<double(void)>>& functions )
{
  auto cov_matrix = covarianceMatrix( functions );
  std::vector<double> diag(functions.size());
  for( size_t i = 0; i < functions.size(); ++i) diag[i] = cov_matrix(i,i);
  for( size_t i = 0; i < functions.size(); ++i ){
    for( size_t j = 0; j < functions.size(); ++j ){
      cov_matrix(i,j) = cov_matrix(i,j) / sqrt( diag[i] * diag[j] );
    }
  }
  return cov_matrix;
}

const std::vector<MinuitParameter*>& LinearErrorPropagator::params() const { return m_parameters; }
