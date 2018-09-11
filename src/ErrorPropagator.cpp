#include "AmpGen/ErrorPropagator.h"

using namespace AmpGen;

GaussErrorPropagator::GaussErrorPropagator( const std::vector<AmpGen::MinuitParameter*>& params, const TMatrixD& reducedCovariance,
    TRandom3* rnd )
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
  transpose();
}

void GaussErrorPropagator::perturb()
{
  const unsigned int N = m_decomposedCholesky.GetNrows();
  TVectorD e( N );
  for ( unsigned int i = 0; i < N; ++i ) e[i] = m_rand->Gaus( 0, 1 );

  TVectorD p = m_decomposedCholesky * e; // perturbation in the usual parameter basis
  for ( int j = 0; j < p.GetNrows(); ++j ) {
    auto f = m_parameters[j];
    f->setCurrentFitVal( m_startingValues[j] + p[j] );
  }
}

void GaussErrorPropagator::reset()
{
  for ( unsigned int j = 0; j < m_parameters.size(); ++j ) m_parameters[j]->setCurrentFitVal( m_startingValues[j] );
}

void GaussErrorPropagator::transpose()
{
  for ( int i = 0; i < m_decomposedCholesky.GetNrows(); ++i ) {
    for ( int j = i + 1; j < m_decomposedCholesky.GetNrows(); ++j ) {
      double tmp = m_decomposedCholesky( j, i );
      m_decomposedCholesky( j, i ) = m_decomposedCholesky( i, j );
      m_decomposedCholesky( i, j ) = tmp;
    }
  }
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

void LinearErrorPropagator::reset()
{
  for ( int x = 0; x < m_cov.GetNrows(); ++x ) {
    for ( int y = 0; y < m_cov.GetNrows(); ++y ) {
      m_cov( x, y ) = 0;
    }
  }
  m_parameters.clear();
}
