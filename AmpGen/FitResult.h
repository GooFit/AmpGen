#ifndef AMPGEN_FITRESULT_H
#define AMPGEN_FITRESULT_H

#include "TMatrixD.h"

#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"

namespace AmpGen
{
  class Minimiser;

  class FitResult
  {
  private:
    double m_chi2;
    double m_LL;
    double m_nBins;
    double m_nParam;
    int m_status;
    EventType m_eventType;
    std::map<std::string, double> m_observables;
    std::vector<FitFraction> m_fitFractions;
    TMatrixD m_covarianceMatrix;
    std::shared_ptr<MinuitParameterSet> m_mps;
    bool m_fitted;
    std::map<std::string, unsigned int> m_covMapping;
    std::string getLastLine( std::ifstream& in ) const;
    void addToParameters( const std::string& line );
    void addToObservables( const std::string& line );
    void setFitQuality( const std::string& line );

  public:
    FitResult( const FitResult& other );

    void addObservable( const std::string& name, const double& F );
    FitResult( const Minimiser& mini );
    FitResult( const MinuitParameterSet& mps, const TMatrixD& covMini );

    FitResult( const std::string& filename, const EventType& evtType = EventType() );
    FitResult() : m_chi2( 0 ), m_LL( -999 ), m_nBins( 0 ), m_nParam( 0 ), m_status( -1 ), m_fitted( 0 ) {}

    bool readFile( const std::string& fname );
    void writeToFile( const std::string& fname );

    double chi2() const { return m_chi2; }
    double LL() const { return m_LL; }
    int status() const { return m_status; }
    int nParam() const { return m_nParam; }
    int nBins() const { return m_nBins; }

    EventType eventType() const { return m_eventType; }
    std::map<std::string, double> observables() const { return m_observables; }
    std::vector<FitFraction> fitFractions() const { return m_fitFractions; }
    std::shared_ptr<MinuitParameterSet> mps() const { return m_mps; }

    double dof() const { return m_nBins - m_nParam - 1; }
    std::vector<MinuitParameter*> getParameters() const;
    std::vector<FitFraction> getFitFractions() const { return m_fitFractions; }
    MinuitParameterSet* MPS() const { return &( *m_mps ); }
    TMatrixD cov() const { return m_covarianceMatrix; }
    double cov( const size_t& x, const size_t& y ) const { return m_covarianceMatrix( x, y ); }

    void print() const;

    std::vector<MinuitParameter*> getFloating( const bool& extended = false ) const;
    TMatrixD getReducedCovariance( const bool& extended = false ) const;
    LinearErrorPropagator getErrorPropagator( const bool& extended = false ) const;
    void addChi2( const double& chi2, const double& nBins );
    void addFractions( const std::vector<FitFraction>& fractions );
    void addFraction( const std::string& name, const double& frac, const double& err )
    {
      m_fitFractions.emplace_back( name, frac, err );
    }
    void clearFitFractions() { m_fitFractions.clear(); }
    void setCov( const size_t& x, const size_t& y, const double& F ) { m_covarianceMatrix( x, y ) = F; }
    double correlation( const std::string& x, const std::string& y )
    {
      auto ix = m_covMapping[x];
      auto iy = m_covMapping[y];
      return m_covarianceMatrix( ix, iy ) / sqrt( m_covarianceMatrix( ix, ix ) * m_covarianceMatrix( iy, iy ) );
    }
  };
} // namespace AmpGen

#endif
