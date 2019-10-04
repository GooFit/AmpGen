#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixTUtils.h>
#include <memory.h>
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

FitResult::FitResult() = default; 

FitResult::FitResult( const FitResult& other )
  : m_mps( std::make_shared<MinuitParameterSet>( *other.mps() ) )
  , m_chi2( other.chi2() )
  , m_LL( other.LL() )
  , m_nBins( other.nBins() )
  , m_nParam( other.nParam() )
  , m_status( other.status() )
  , m_observables( other.observables() )
  , m_fitFractions( other.fitFractions() )
  , m_covarianceMatrix(other.cov())
{
}

FitResult::FitResult( const std::string& filename ) :
  m_mps( std::make_shared<MinuitParameterSet>() ) 
{
  m_fitted = readFile( filename );
}

FitResult::FitResult( const Minimiser& mini )
  : m_mps  ( std::make_shared<MinuitParameterSet>( *mini.parSet() ) )
  , m_LL   ( mini.FCN() )
  , m_nParam( 0 )
  , m_status( mini.status() )
  , m_covarianceMatrix( mini.covMatrixFull() ) {
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->at(i)->isFree() ) m_nParam++;
    m_covMapping[ m_mps->at(i)->name() ] = i;
  }
}

FitResult::FitResult( const MinuitParameterSet& mps, const TMatrixD& covMini ) :
  m_mps( std::make_shared<MinuitParameterSet>( mps ) )
{
  if ( int( mps.size() ) != covMini.GetNcols() ) {
    ERROR( "Minuit parameter set size does not match covariance matrix size!" );
  }
  m_covarianceMatrix.ResizeTo( covMini.GetNcols(), covMini.GetNrows() );
  m_covarianceMatrix = covMini;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->at(i)->isFree()) m_nParam++;
  }
}

bool FitResult::readFile( const std::string& fname )
{
  std::ifstream checkIsClosed( fname );
  if ( !checkIsClosed.is_open() 
      || checkIsClosed.peek() == std::ifstream::traits_type::eof() ) return false;
  checkIsClosed.close();
  auto lines = vectorFromFile(fname);
  if( *lines.rbegin() != "End Log" ){
    ERROR("File not properly closed: " << *lines.rbegin() );
    return false; 
  }
  std::vector<std::string> parameterLines;
  for( auto& line : lines ){
    const std::string name = split( line, ' ' )[0];
    if ( name == "Parameter" ) parameterLines.push_back( line );
    else if ( name == "FitQuality" )  this->setFitQuality( line );
    else if ( name == "FitFraction" ) this->m_fitFractions.emplace_back(line);
    else if ( name == "Observable" )  this->addToObservables( line );
  }
  size_t nParameters = parameterLines.size();
  m_covarianceMatrix.ResizeTo( parameterLines.size(), parameterLines.size() );
  for (size_t i = 0; i < nParameters; ++i ) {
    auto tokens             = split( parameterLines[i], ' ' );
    m_covMapping[tokens[1]] = i;
    m_mps->add( new MinuitParameter( tokens[1], parse<Flag>(tokens[2]), stod( tokens[3] ), stod( tokens[4] ), 0, 0 ) );
    for (size_t j = 0; j < nParameters; ++j ) m_covarianceMatrix( i, j ) = stod( tokens[5 + j] );
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
  for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
    auto param = m_mps->at(i);
    outlog << "Parameter"
      << " " << param->name() << " " << to_string<Flag>(param->flag()) << " " << param->mean() << " "
      << ( param->isFree() ? m_mps->at(i)->err() : 0 ) << " ";
    for (size_t j = 0; j < (size_t)m_covarianceMatrix.GetNcols(); ++j ) outlog << m_covarianceMatrix[i][j] << " ";
    outlog << std::endl;
  }
  outlog << std::setprecision( 8 );
  outlog << "FitQuality " << m_chi2 << " " << m_nBins << " " << m_nParam << " " << m_LL    << " " << m_status << "\n";
  for ( auto& f : m_fitFractions )  outlog << "FitFraction " << f.name() << " " << f.val() << " " << f.err()  << "\n";
  for ( auto& o : m_observables )   outlog << "Observable "  << o.first  << " " << o.second << "\n";
  outlog << "End Log\n";
  outlog.close();
}

void FitResult::print() const
{
  INFO( "Chi2 per bin = " << m_chi2 / m_nBins );
  INFO( "Chi2 per dof = " << m_chi2 / dof() );
  INFO( "-2LL         = " << m_LL );
  INFO( "Fit Status   = " << m_status );
}

std::vector<MinuitParameter*> FitResult::parameters() const
{
  std::vector<MinuitParameter*> params( m_mps->size() );
  std::copy( m_mps->begin(), m_mps->end(), params.begin() );
  return params;
}

std::vector<MinuitParameter*> FitResult::floating( const bool& extended ) const
{
  std::vector<MinuitParameter*> floating;
  for ( auto& param : *m_mps ) {
    if ( ( param->isFree() || extended ) && param->err() > 1e-6 ) floating.push_back( param );
  }
  return floating;
}

TMatrixD FitResult::getReducedCovariance( const bool& extended ) const
{
  std::vector<unsigned int> floating_indices;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    auto param = m_mps->at(i);
    if ( ( param->isFree() || extended ) && param->err() > 1e-6 ) floating_indices.push_back( i );
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
  return LinearErrorPropagator( getReducedCovariance( extended ), floating( extended ) );
}

void FitResult::addChi2( const double& chi2, const double& nBins )
{
  m_chi2  = chi2;
  m_nBins = nBins;
}

void FitResult::addFractions( const std::vector<FitFraction>& fractions ) 
{ 
  m_fitFractions = fractions; 
}

double FitResult::chi2() const { return m_chi2; }
double FitResult::LL() const { return m_LL; }
int FitResult::status() const { return m_status; }
int FitResult::nParam() const { return m_nParam; }
int FitResult::nBins() const { return m_nBins; }

std::map<std::string, double> FitResult::observables() const { return m_observables; }
std::shared_ptr<MinuitParameterSet> FitResult::mps() const { return m_mps; }

double FitResult::dof() const { return m_nBins - m_nParam - 1; }
std::vector<FitFraction> FitResult::fitFractions() const { return m_fitFractions; }
TMatrixD FitResult::cov() const { return m_covarianceMatrix; }
double FitResult::cov( const size_t& x, const size_t& y ) const { return m_covarianceMatrix( x, y ); }
double FitResult::cov( const std::string& x, const std::string& y ) const { return m_covarianceMatrix( m_covMapping.at(x), m_covMapping.at(y) ); }

void FitResult::addFraction( const std::string& name, const double& frac, const double& err )
{
  m_fitFractions.emplace_back( name, frac, err );
}
void FitResult::clearFitFractions() { m_fitFractions.clear(); }
void FitResult::setCov( const size_t& x, const size_t& y, const double& F ) { m_covarianceMatrix( x, y ) = F; }
double FitResult::correlation( const std::string& x, const std::string& y ) const
{
  auto tx = m_covMapping.find(x);
  auto ty = m_covMapping.find(y);
  if( tx == m_covMapping.end() || ty == m_covMapping.end() ){
    ERROR("Parameter not found: " << x << ", " << y );
    return -1;
  }
  return cov(tx->second, ty->second)/sqrt(cov(tx->second, tx->second)*cov(ty->second, ty->second));
}
