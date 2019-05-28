#ifndef AMPGEN_ERRORPROPAGATOR_H
#define AMPGEN_ERRORPROPAGATOR_H

#include <cmath>
#include <stddef.h>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <array>
#include <utility>

#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

#include "TDecompChol.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include "TVectorDfwd.h"
#include "TVectorT.h"

class TRandom3;

namespace AmpGen
{
  class MinuitParameterSet;
  /** @class LinearErrorPropagator
   Propagates uncertainties on functors using either a MinuitParameterSet (thus assuming a diagonal covariance matrix) 
   or with a set of parameters and a covariance matrix (i.e. the product of a fit). 
   Can be used for calculating the uncertainties on observables such as Fit fractions, 
   or for producing uncertainty bands in predictions on figures. 
   Usage is 
  */
  class LinearErrorPropagator
  {
    public:
      ///< Constructor for LinearErrorPropagator taking a vector of free parameters, assumes a diagonal covariance matrix taking the parameter uncertainties.
      LinearErrorPropagator( const std::vector<MinuitParameter*>& params );
      ///< Constructor for LinearErrorPropagator, taking a MinuitParameterSet as argument, assumes a diagonal coviarance matrix using the uncertainties on parameters.
      LinearErrorPropagator( const MinuitParameterSet& params );
      ///< Constructor for LinearErrorPropagator, taking a covariance matrix and a vector parameters
      LinearErrorPropagator( const TMatrixD& reducedCovarianceMatrix, const std::vector<MinuitParameter*>& params );
      
      ///< Calculates the uncertainty on functor fcn.  
      template <class FCN> double operator()( FCN fcn ) const
      {
        return getError( std::function<double(void)>(fcn) );
      }

      ///< Calculates the error on functor fcn (should take no arguments and return a double).
      double getError( const std::function<double(void)>& fcn ) const;

      ///< Calculate the uncertainties on a functor that returns a vector of size RANK. 
      std::vector<double> getVectorError( const std::function<std::vector<double>(void)>& fcn, size_t RANK ) const;

      /// Calculate the variance-covariance matrix of functor set functions. 
      std::vector<double> combinationWeights( const std::vector<std::function<double(void)>>& functions );
      TMatrixD covarianceMatrix(const std::vector<std::function<double(void)>>& functions );

      /// Calculate the correlation matrix of functor set functions. 
      TMatrixD correlationMatrix(const std::vector<std::function<double(void)>>& functions );

      /// Calculate the covariance-matrix weighted combination of functor set functions. 
      std::pair<double,double> combinationCovWeighted(const std::vector<std::function<double(void)>>& functions );

      void add( const LinearErrorPropagator& p2 );
      void reset(); 
      size_t size()                                 const;
      const TMatrixD& cov()                         const;
      const std::map<std::string, size_t> posMap()  const;
      const std::vector<MinuitParameter*>& params() const; 

    private:
      TMatrixD                              m_cov;
      std::vector<MinuitParameter*>         m_parameters;

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
  };

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

} // namespace AmpGen
#endif
