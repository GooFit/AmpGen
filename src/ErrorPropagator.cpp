#include "AmpGen/ErrorPropagator.h"

#include <ostream>

#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Minimiser.h"
#include "TDecompChol.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TVectorT.h"
#include "TFile.h"
#include "TTree.h"

typedef TVectorT<double> TVectorD; 

using namespace AmpGen;

GaussErrorPropagator::GaussErrorPropagator( const TMatrixD& reducedCovariance, const std::vector<MinuitParameter*>& params, TRandom3* rnd )
  : m_parameters( params ), m_rand( rnd ), m_decomposedCholesky( params.size(), params.size() )
{
  for ( size_t x = 0; x < params.size(); ++x ) {
    auto p = params[x];
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

LinearErrorPropagator::LinearErrorPropagator( Minimiser* mini )
  : m_cov( mini->covMatrix() ) 
{
  for( auto& param : *mini->parSet() ){
    if( !param->isFree() ) continue;
    m_parameters.push_back( param );
  }
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
    // CHECK_BLINDING
    DEBUG( "Perturbing parameter: [" << m_parameters[i]->name() << "] " << m_parameters[i]->mean() << " by "
        << sqrt( m_cov( i, i ) ) << " " << m_parameters[i] );
    errorVec(i) = derivative(fcn,i);
    fcn(); 
  }
  return sqrt( Dot( errorVec, m_cov * errorVec ) );
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
    // CHECK_BLINDING
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
  for ( unsigned int j = 0; j < RANK; ++j ) rt[j] = sqrt( Dot( errorVec[j] , m_cov * errorVec[j] ) );
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

NonlinearErrorPropagator::NonlinearErrorPropagator( Minimiser* mini ) : m_mini(mini) {}

template <typename target_t, typename proposal_t> 
class MetropolisHastings {
  public:

    MetropolisHastings( const target_t& target, const proposal_t& proposal, const TVectorD& errors )
      : target(target), proposal(proposal),errors(errors) {}
  
  TVectorD operator() () const
  {
    auto x = proposal(); 
    for( int i = 0 ; i != 500; ++i )
    {
      auto N = x.GetNrows();
      auto x_prime = TVectorD( x.GetNrows() ); 
      for( int i = 0 ; i != N; ++i ) x_prime(i) = x(i) + errors(i) * gRandom->Gaus(); 
      auto r = target( x_prime ) /target(x);
      x = r > gRandom->Uniform() ? x_prime : x; 
    }
    return x; 
  }

  target_t target;
  proposal_t proposal;
  TVectorD errors;
};

template <typename target_t, typename proposal_t> MetropolisHastings<target_t, proposal_t> make_sampler ( const target_t& target, const proposal_t& proposal, const TVectorD& err)
{
  return MetropolisHastings<target_t, proposal_t> ( target, proposal, err);
}

TMatrixD NonlinearErrorPropagator::correlationMatrix( const std::vector<std::function<double(void)>>& functions, TRandom3* rnd, const unsigned& nSamples) 
{
  std::vector<std::pair<MinuitParameter*, double>> parameters;
  std::vector<double> function_averages; 
  for( auto& function : functions ) function_averages.push_back(  function() );

  for ( auto& param : *m_mini->parSet() )
  {
    if( param->isFree() ) parameters.emplace_back(param, param->mean());  
  }
  double L0 = m_mini->FCN(); 

  auto reducedCovariance = m_mini->covMatrix(); 
  TMatrixD invCovariance     = reducedCovariance;
  invCovariance.Invert(); 
  TDecompChol decomposed( reducedCovariance );
  decomposed.Decompose();
  auto decomposedCholesky = decomposed.GetU();
  for ( int i = 0; i < decomposedCholesky.GetNrows(); ++i ) {
    for ( int j = i + 1; j < decomposedCholesky.GetNrows(); ++j ){
      std::swap( decomposedCholesky(i, j), decomposedCholesky(j, i)  );
    }
  }
  auto Q = []( const auto& x, const auto& A ){ return ( x * A * TMatrixD( TMatrixD::kTransposed , x ) )(0,0); };
  std::vector<double> f_prime(functions.size());
  
  std::vector<double> parameters_v( parameters.size() );
  const unsigned int N = decomposedCholesky.GetNrows();  
  double weight = 1 ;
  TTree* tree = nullptr; 
  TFile* file = TFile::Open("test.root","RECREATE");
  tree = new TTree("t","t");
  
  for( size_t i = 0 ; i != functions.size(); ++i ) tree->Branch( ("f_"+std::to_string(i)).c_str(),  &f_prime[i] ); 
  
  for( unsigned int i = 0 ; i != N; ++i ) tree->Branch( parameters[i].first->name().c_str(), &parameters_v[i] ); 
  TVectorD p0(N);
  TVectorD err(N);
  for(unsigned i = 0 ; i != N; ++i ){
    p0(i) = parameters[i].second; 
    err(i) = parameters[i].first->err(); 
  }
  auto proposal_distribution = [&]()
  {
    TVectorD e( N );
    for ( unsigned i = 0; i != N; ++i ) e[i] = rnd->Gaus( 0, 1 );
    return p0 + decomposedCholesky * e;
  };
  auto target_distribution = [&]( TVectorD& state )
  {
    for( unsigned int i = 0 ; i!= N; ++i ) parameters[i].first->setCurrentFitVal( state(i ) );
    return exp( -0.5*(m_mini->FCN() -L0 ) );
  };

  TMatrixD rt( functions.size(), functions.size() );
  double w_sum = 0; 
  auto mh = make_sampler(target_distribution, proposal_distribution, err);    
  for( unsigned sample = 0 ; sample != nSamples ; ++sample )
  {
    if( sample % 100 == 0 ) INFO( sample << " / " << nSamples << " completed"); 
    auto p = mh(); 
    for( unsigned int i = 0 ; i!= N; ++i ) parameters[i].first->setCurrentFitVal( p(i ) );
    for( unsigned int i = 0 ; i!= N; ++i ) parameters_v[i] = p(i);
     
    for( unsigned i = 0 ; i != functions.size(); ++i ) f_prime[i] = ( functions[i]() ); 

    for( unsigned i = 0 ; i != functions.size(); ++i ) 
    {
      for( unsigned j = 0 ; j != functions.size(); ++j )
      {
        weight = 1; 
        rt(i,j) += weight * ( f_prime[i] - function_averages[i] ) * ( f_prime[j] - function_averages[j] ); 
        w_sum += weight; 
      }
    }
    if( tree != nullptr ) tree->Fill();
  }
  tree->Write();
  file->Close();
  rt *= 1./w_sum; 
  return rt; 
}


