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
    GaussErrorPropagator( const TMatrixD& reducedCovariance,
                          const std::vector<MinuitParameter*>& params, 
                          TRandom3* rnd );
    void perturb();
    void reset();
    void transpose();

  private:
    std::vector<MinuitParameter*>         m_parameters;
    std::vector<double>                   m_startingValues;
    TRandom3*                             m_rand;
    TMatrixD                              m_decomposedCholesky;
  };

  class LinearErrorPropagator
  {

  public:
    LinearErrorPropagator( const std::vector<MinuitParameter*>& params );
    LinearErrorPropagator( const MinuitParameterSet& params );

    LinearErrorPropagator( const TMatrixD& reducedCovarianceMatrix,
                           const std::vector<MinuitParameter*>& params );
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

    template <class FCN> double getError( FCN fcn ) const
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
    template <size_t RANK, class FCN> std::array<double, RANK> getVectorError( FCN fcn ) const
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
    template <class FCN> std::vector<double> getVectorError( FCN fcn, size_t RANK ) const
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

    TMatrixD propagatedCovarianceMatrix(const std::vector<std::function<double(void)>>& functions ) ;
    std::pair<double,double> combinationCovWeighted(const std::vector<std::function<double(void)>>& functions ) ;

    template <class FCN> double operator()( FCN fcn ) const
    {
      return getError( fcn );
    }

    void add( const LinearErrorPropagator& p2 );

    void reset();
   
    size_t size() const ;
    const TMatrixD& cov() const ;
    const std::map<std::string, size_t> posMap() const;
    const std::vector<MinuitParameter*>& params() const; 

  private:
    TMatrixD                              m_cov;
    std::vector<MinuitParameter*>         m_parameters;
  };
} // namespace AmpGen
#endif
