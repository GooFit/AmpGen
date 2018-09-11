#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixTUtils.h>
#include <iomanip>
#include <istream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

FitResult::FitResult( const std::string& filename, const EventType& evtType )
{
  m_eventType = evtType;
  m_fitted    = readFile( filename );
}

std::string FitResult::getLastLine( std::ifstream& in ) const
{
  std::string line;
  while ( in >> std::ws && std::getline( in, line ) )
    ;
  return line;
}

bool FitResult::readFile( const std::string& fname )
{
  std::ifstream CHECK( fname );
  if ( !CHECK.is_open() || CHECK.peek() == std::ifstream::traits_type::eof() ) {
    return false;
  }
  if ( getLastLine( CHECK ) != "End Log" ) {
    ERROR( "File not properly close " << fname );
    return false;
  }
  CHECK.close();

  // auto lines = vectorFromFile( fname );
  std::vector<std::string> parameterLines;

  processFile( fname, [this, &parameterLines]( auto& line ) {
    const std::string name = split( line, ' ' )[0];
    if ( name == "Parameter" )
      parameterLines.push_back( line );
    else if ( name == "FitQuality" )
      this->setFitQuality( line );
    else if ( name == "FitFraction" )
      this->m_fitFractions.emplace_back( line, m_eventType );
    else if ( name == "Observable" )
      this->addToObservables( line );
  } );

  unsigned int nParameters = parameterLines.size();

  m_covarianceMatrix.ResizeTo( parameterLines.size(), parameterLines.size() );
  m_mps = std::make_shared<MinuitParameterSet>();
  for ( unsigned int i = 0; i < nParameters; ++i ) {
    auto tokens             = split( parameterLines[i], ' ' );
    m_covMapping[tokens[1]] = i;
    m_mps->add( new MinuitParameter( tokens[1], MinuitParameter::Flag(stoi( tokens[2] ) ), stod( tokens[3] ), stod( tokens[4] ), 0, 0 ) );
    for ( unsigned int j = 0; j < nParameters; ++j ) m_covarianceMatrix( i, j ) = stod( tokens[5 + j] );
  }

  return true;
}

void FitResult::setFitQuality( const std::string& line )
{
  auto tokens = split( line, ' ' );
  bool status=true;
  if( tokens.size() != 6 ){
   WARNING("Cannot pass FitQuality line: " << line );
   return;
  }
  m_chi2      = lexical_cast<double>( tokens[1] , status );
  m_nBins     = lexical_cast<double>( tokens[2] , status );
  m_nParam    = lexical_cast<double>( tokens[3] , status );
  m_LL        = lexical_cast<double>( tokens[4] , status );
  m_status    = lexical_cast<int>   ( tokens[5] , status );
}
void FitResult::addToObservables( const std::string& line )
{
  auto tokens              = split( line, ' ' );
  m_observables[tokens[1]] = stod( tokens[2] );
}

void FitResult::addObservable( const std::string& name, const double& F ) { m_observables[name] = F; }

void FitResult::writeToFile( const std::string& fname )
{
  std::ofstream outlog;
  outlog.open( fname );

  for ( int i = 0; i < m_covarianceMatrix.GetNrows(); ++i ) {
    auto param = m_mps->getParPtr( i );
    outlog << "Parameter"
           << " " << param->name() << " " << param->iFixInit() << " " << param->mean() << " "
           << m_mps->getParPtr( i )->err() << " ";
    for ( int j = 0; j < m_covarianceMatrix.GetNcols(); ++j ) outlog << m_covarianceMatrix[i][j] << " ";
    outlog << std::endl;
  }
  outlog << std::setprecision( 8 );
  outlog << "FitQuality " << m_chi2 << " " << m_nBins << " " << m_nParam << " " << m_LL << " " << m_status << "\n";
  for ( auto& p : m_fitFractions ) {
    outlog << "FitFraction " << p.name() << " " << p.val() << " " << p.err() << "\n";
  }
  for ( auto& ob : m_observables ) {
    outlog << "Observable " << ob.first << " " << ob.second << "\n";
  }
  outlog << "End Log\n";
  outlog.close();
}

FitResult::FitResult( const AmpGen::Minimiser& mini )
{
  m_mps  = std::make_shared<AmpGen::MinuitParameterSet>( *mini.parSet() );
  m_LL   = mini.FCN();
  auto M = mini.covMatrixFull();
  m_covarianceMatrix.ResizeTo( M.GetNcols(), M.GetNrows() );
  m_covarianceMatrix = M;
  m_status           = mini.status();
  m_nParam           = 0;
  for ( unsigned int i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->getParPtr( i )->iFixInit() == 0 ) m_nParam++;
    m_covMapping[ m_mps->getParPtr(i)->name() ] = i;
  }
}

FitResult::FitResult( const AmpGen::MinuitParameterSet& mps, const TMatrixD& covMini ) : FitResult()
{
  m_mps = std::make_shared<AmpGen::MinuitParameterSet>( mps );
  if ( int( mps.size() ) != covMini.GetNcols() ) {
    ERROR( "Minuit parameter set size does not match covariance matrix size!" );
  }
  m_covarianceMatrix.ResizeTo( covMini.GetNcols(), covMini.GetNrows() );
  m_covarianceMatrix = covMini;
  for ( unsigned int i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->getParPtr( i )->iFixInit() == 0 ) m_nParam++;
  }
}

void FitResult::print() const
{
  INFO( "Chi2 per bin = " << m_chi2 / m_nBins );
  INFO( "Chi2 per dof = " << m_chi2 / dof() );
  INFO( "-2LL         = " << m_LL );
  INFO( "Fit Status   = " << m_status );
}

FitResult::FitResult( const FitResult& other )
    : m_chi2( other.chi2() )
    , m_LL( other.LL() )
    , m_nBins( other.nBins() )
    , m_nParam( other.nParam() )
    , m_status( other.status() )
    , m_eventType( other.eventType() )
    , m_observables( other.observables() )
    , m_fitFractions( other.fitFractions() )
    , m_mps( std::make_shared<MinuitParameterSet>( *other.mps() ) )
{
  m_covarianceMatrix.ResizeTo( other.cov().GetNrows(), other.cov().GetNcols() );
  m_covarianceMatrix = other.cov();
}

std::vector<AmpGen::MinuitParameter*> FitResult::getParameters() const
{
  std::vector<AmpGen::MinuitParameter*> params;
  for ( auto& param : *m_mps ) params.push_back( param );
  return params;
}

std::vector<AmpGen::MinuitParameter*> FitResult::getFloating( const bool& extended ) const
{
  std::vector<AmpGen::MinuitParameter*> floating;
  for ( auto& param : *m_mps ) {
    if ( ( param->iFixInit() == 0 || extended ) && param->err() > 1e-6 ) floating.push_back( param );
  }
  /*
  if( extended ){
    DEBUG("Got extended error propagator:");
    for( auto& param : floating ){
      INFO( param->name() );
    }
  }
  */
  return floating;
}

TMatrixD FitResult::getReducedCovariance( const bool& extended ) const
{
  std::vector<unsigned int> floating_indices;
  for ( unsigned int i = 0; i < m_mps->size(); ++i ) {
    auto param = m_mps->getParPtr( i );
    if ( ( param->iFixInit() == 0 || extended ) && param->err() > 1e-6 ) floating_indices.push_back( i );
  }
  TMatrixD reducedCov( floating_indices.size(), floating_indices.size() );
  if ( int( floating_indices.size() ) > m_covarianceMatrix.GetNrows() ) {
    ERROR( "Looking for more floating indices than real indices: " << m_covarianceMatrix.GetNrows() << " "
                                                                   << floating_indices.size() );
    return reducedCov;
  }
  for ( unsigned int i = 0; i < floating_indices.size(); ++i ) {
    for ( unsigned int j = 0; j < floating_indices.size(); ++j ) {
      reducedCov( i, j ) = m_covarianceMatrix( floating_indices[i], floating_indices[j] );
    }
  }
  return reducedCov;
}

LinearErrorPropagator FitResult::getErrorPropagator( const bool& extended ) const
{
  return LinearErrorPropagator( getReducedCovariance( extended ), getFloating( extended ) );
}

void FitResult::addChi2( const double& chi2, const double& nBins )
{
  m_chi2  = chi2;
  m_nBins = nBins;
}

void FitResult::addFractions( const std::vector<FitFraction>& fractions ) { m_fitFractions = fractions; }
