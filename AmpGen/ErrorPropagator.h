#ifndef AMPGEN_ERRORPROPAGATOR_H
#define AMPGEN_ERRORPROPAGATOR_H

#include <string>
#include <vector>
#include <map>
#include <functional>

#include "TDecompChol.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TVectorD.h"

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

namespace AmpGen
{
  class GaussErrorPropagator
  {
  public:
    GaussErrorPropagator( const std::vector<AmpGen::MinuitParameter*>& params, const TMatrixD& reducedCovariance,
                          TRandom3* rnd );
    void perturb();
    void reset();
    void transpose();

  private:
    std::vector<AmpGen::MinuitParameter*> m_parameters;
    std::vector<double> m_startingValues;
    TRandom3* m_rand;
    TMatrixD m_decomposedCholesky;
  };

  class LinearErrorPropagator
  {

  public:
    LinearErrorPropagator( const TMatrixD& reducedCovarianceMatrix,
                           const std::vector<AmpGen::MinuitParameter*>& params )
        : m_cov( reducedCovarianceMatrix ), m_parameters( params )
    {
    }
    template <class FCN> double derivative( FCN fcn, const size_t& i ) const
    {
      double startingValue = m_parameters[i]->mean();
      m_parameters[i]->setCurrentFitVal( startingValue + sqrt( m_cov( i, i ) ) );
      double plus_variation = fcn();
      m_parameters[i]->setCurrentFitVal( startingValue - sqrt( m_cov( i, i ) ) );
      double minus_variation = fcn();
      m_parameters[i]->setCurrentFitVal( startingValue );
      return ( plus_variation - minus_variation ) / ( 2 * sqrt( m_cov( i, i ) ) );
    }

    template <class FCN>
    double getError( FCN fcn ) const
    {
      unsigned int N = m_cov.GetNrows();
      TVectorD errorVec( N );
      for ( unsigned int i = 0; i < N; ++i ) {
        DEBUG( "Perturbing parameter: [" << m_parameters[i]->name() << "] " << startingValue << " by "
                                         << sqrt( m_cov( i, i ) ) << " " << m_parameters[i] );
        errorVec(i) = derivative(fcn,i);
        fcn(); /// ensure any cache state is okay still ///
      }
      return sqrt( errorVec * ( m_cov * errorVec ) );
    }
    template <size_t RANK, class FCN>
    std::array<double, RANK> getVectorError( FCN fcn ) const
    {
      unsigned int N = m_cov.GetNrows();
      std::array<TVectorD, RANK> errorVec;
      for ( unsigned int i = 0; i < RANK; ++i ) errorVec[i].ResizeTo( N );
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
      std::array<double, RANK> rt;
      for ( unsigned int j = 0; j < RANK; ++j ) rt[j] = sqrt( errorVec[j] * ( m_cov * errorVec[j] ) );

      return rt;
    }
    template <class FCN>
    std::vector<double> getVectorError( FCN fcn, size_t RANK ) const
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

    TMatrixD propagatedCovarianceMatrix(const std::vector<std::function<double(void)>>& functions ){
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


    template <class FCN>
    double operator()( FCN fcn ) const
    {
      return getError( fcn );
    }

    void add( const LinearErrorPropagator& p2 );

    void reset();
   
    size_t size() const { return m_parameters.size(); }
    const TMatrixD& cov() const { return m_cov; }
    const std::map<std::string, size_t> posMap() const
    {
      std::map<std::string, size_t> pMap;
      for ( size_t pos = 0; pos != m_parameters.size(); ++pos ) {
        pMap[m_parameters[pos]->name()] = pos;
      }
      return pMap;
    }
    const std::vector<AmpGen::MinuitParameter*>& params() const { return m_parameters; }

  private:
    TMatrixD m_cov;
    std::vector<AmpGen::MinuitParameter*> m_parameters;
  };
} // namespace AmpGen
#endif
