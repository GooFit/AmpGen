#ifndef AMPGEN_FITRESULT_H
#define AMPGEN_FITRESULT_H

#include "TMatrixD.h"

#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Utilities.h"

namespace AmpGen
{
  class Minimiser;
  class LinearErrorPropagator; 

  class FitResult
  {
  public:
    FitResult(); 
    FitResult( const FitResult& other );
    FitResult( const std::string& filename );
    FitResult( const Minimiser& mini );
    FitResult( const MinuitParameterSet& mps, const TMatrixD& covMini );

    void addObservable( const std::string& name, const double& F );
    void addChi2( const double& chi2, const double& nBins );
    void addFractions( const std::vector<FitFraction>& fractions );
    void addFraction( const std::string& name, const double& frac, const double& err );
    void setCov( const size_t& x, const size_t& y, const double& F );
    void writeToFile( const std::string& fname );
    void clearFitFractions();

    bool readFile( const std::string& fname );

    double chi2() const;
    double LL()   const;
    double dof()  const;
    double cov(const size_t& x, const size_t& y ) const;
    double cov(const std::string& x, const std::string& y ) const; 
    double correlation( const std::string& x, const std::string& y ) const;
    
    int status() const;
    int nParam() const;
    int nBins()  const; 

    std::map<std::string, double> observables() const;
    std::shared_ptr<MinuitParameterSet> mps()   const;

    std::vector<FitFraction> fitFractions()     const;
    std::vector<MinuitParameter*> parameters()  const;
    std::vector<MinuitParameter*> floating(const bool& extended = false) const;
    TMatrixD cov() const;
    
    void print() const;

    TMatrixD getReducedCovariance( const bool& extended = false ) const;
    LinearErrorPropagator getErrorPropagator( const bool& extended = false ) const;

  private:
    std::shared_ptr<MinuitParameterSet> m_mps;
    double                              m_chi2   = {0};
    double                              m_LL     = {-999};
    double                              m_nBins  = {0};
    double                              m_nParam = {0};
    int                                 m_status = {-1};
    bool                                m_fitted = {false};
    std::map<std::string, double>       m_observables;
    std::vector<FitFraction>            m_fitFractions;
    TMatrixD                            m_covarianceMatrix;
    std::map<std::string, unsigned int> m_covMapping;

    void addToParameters( const std::string& line );
    void addToObservables( const std::string& line );
    void setFitQuality( const std::string& line );
  };
} // namespace AmpGen

#endif
