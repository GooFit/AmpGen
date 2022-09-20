#include <memory.h>
#include <iomanip>
#include <istream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <set>

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
  : m_mps( other.mps()  )
  , m_chi2( other.chi2() )
  , m_LL( other.LL() )
  , m_Edm( other.Edm() )
  , m_NCalls( other.NCalls() )
  , m_nBins( other.nBins() )
  , m_nParam( other.nParam() )
  , m_status( other.status() )
  , m_observables( other.observables() )
  , m_fitFractions( other.fitFractions() )
  , m_covarianceMatrix(other.cov())
{
}

FitResult::FitResult( const std::string& filename ) :
  m_mps( new MinuitParameterSet() ) 
{
  m_fitted = readFile( filename );
}

FitResult::FitResult( const Minimiser& mini )
  : m_mps  ( mini.parSet() )
  , m_LL   ( mini.FCN() )
  , m_Edm  ( mini.Edm() )
  , m_NCalls ( mini.NCalls() )
  , m_nParam( 0 )
  , m_status( mini.status() )
  , m_covarianceMatrix( mini.covMatrixFull() ) {
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->at(i)->isFree() ||  m_mps->at(i)->isBlind()  ) m_nParam++;
    m_covMapping[ m_mps->at(i)->name() ] = i;
  }
}

FitResult::FitResult( MinuitParameterSet& mps, const TMatrixD& covMini ) : m_mps(&mps)
{
  if ( int( mps.size() ) != covMini.GetNcols() ) {
    ERROR( "Minuit parameter set size does not match covariance matrix size!" );
  }
  m_covarianceMatrix.ResizeTo( covMini.GetNcols(), covMini.GetNrows() );
  m_covarianceMatrix = covMini;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    m_covMapping[m_mps->at(i)->name()] = i; 
    if ( m_mps->at(i)->isFree() || m_mps->at(i)->isBlind() ) m_nParam++;
  }
}

bool FitResult::readFile( const std::string& fname )
{
  auto to_vector_of_doubles = [](auto begin, auto end){ 
    std::vector<double> rt; 
    std::transform(begin, end, std::back_inserter(rt), [](const auto& str){ return stod(str); } ); 
    return rt; };
    auto to_parameter = [&](auto& tokens) -> MinuitParameter* {
    if( tokens.size() != 9 )
    {
      WARNING("Unrecognised number of tokens: " << tokens.size() );
      return nullptr;
    }
    auto args = to_vector_of_doubles( tokens.begin() + 3, tokens.end() ); 
    auto p = new MinuitParameter( tokens[1], parse<Flag>(tokens[2]), args[0], args[1], args[5], args[6] );
    p->setResult(args[0], args[1], args[2], args[3]);
    return p;
  }; 
  std::vector<std::vector<double>> covarianceLines; 
  bool isInCovariance = false; 
  unsigned j = 0; 
  processFile( fname, [&isInCovariance, &j, this, to_parameter, to_vector_of_doubles, &covarianceLines](auto& line ) mutable {
    const auto tokens = split( line, ' ' );
    const auto name = tokens[0]; 
    if ( name == "Parameter" ){ 
      auto param = to_parameter(tokens); 
      if( param != nullptr ) { this->m_mps->add( param ); this->m_covMapping[tokens[1]] = j++; }
    } 
    else if ( name == "FitQuality"  )  this->setFitQuality( line );
    else if ( name == "FitFraction" )  this->m_fitFractions.emplace_back(line);
    else if ( name == "Observable"  )  this->addToObservables( line );
    else if ( name == "ParamsOrder" ){  isInCovariance = false; } 
    else if ( name == "CovarianceMatrix"  ){  isInCovariance = true; return; }
    if( isInCovariance ){
      auto v = to_vector_of_doubles(tokens.begin(), tokens.end() ); 
      covarianceLines.push_back( to_vector_of_doubles(tokens.begin(), tokens.end() ) ); 
    }
  } ); 
  m_covarianceMatrix.ResizeTo(j,j); 
  for (size_t i = 0; i < m_covarianceMatrix.GetNcols(); ++i ) {
    for (size_t j = 0; j < m_covarianceMatrix.GetNrows(); ++j ) m_covarianceMatrix(i, j) = covarianceLines[i][j];
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
  outlog << std::setprecision( 8 );
//   outlog<<"Parameter Name, Flag, mean, error, errNeg, errPos\n";
  for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
    auto param = m_mps->at(i);
    if (param->name().find("_blind") != std::string::npos ) continue;
    double offset = 0; 
    if ( param->isBlind() ){
      offset  =  m_mps->at(param->name()+"_blind")->mean();
      if (!offset ){
	WARNING("Attempting to print the unblind value of a blind parameter. Will skip the parameter.");
	continue;
      }
    }
    outlog << mysprintf("Parameter %-50s %-7s %-12lf %-12lf %-12lf %-12lf %-12lf %-12lf \n", 
        param->name().c_str(), to_string<Flag>(param->flag()).c_str(), 
        param->mean() + offset, param->err(), 
        param->errNeg(), param->errPos(), param->minInit(), param->maxInit() ); 
  }

  outlog << "FitQuality " << m_chi2 << " " << m_nBins << " " << m_nParam << " " << m_LL    << " " << m_status << "\n";
  for ( auto& f : m_fitFractions )  outlog << "FitFraction " << f.name() << " " << f.val() << " " << f.err()  << "\n";
  for ( auto& o : m_observables )   outlog << "Observable "  << o.first  << " " << o.second << "\n";

  outlog << "\nCovarianceMatrix\n";
  for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
    for (size_t j = 0; j < (size_t)m_covarianceMatrix.GetNcols(); ++j ) outlog << m_covarianceMatrix[i][j] << " ";
    outlog<<"\n";
  }
  outlog<<"ParamsOrder ";
  for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
    auto param_i = m_mps->at(i);
    outlog<<param_i->name()<<" "<<to_string<Flag>(param_i->flag())<<" " ;
  }
  outlog << "\nEnd Log\n";
  outlog.close();
}

void FitResult::print() const
{
  INFO( "Chi2 per bin = " << m_chi2 / m_nBins );
  INFO( "Chi2 per dof = " << m_chi2 / dof() );
  INFO( "-2LL         = " << m_LL );
  INFO( "Fit Status   = " << m_status );
  INFO( "NCalls       = " << m_NCalls );
  INFO( "Edm          = " << m_Edm );
  /*
  std::cout<<"\n"<<std::endl;
 
  unsigned longest_parameter_name = 10;
  for( unsigned i = 0 ; i != (unsigned)m_covarianceMatrix.GetNrows(); ++i )
  {
    auto param = m_mps->at(i); 
    if( param->name().size() > longest_parameter_name ) longest_parameter_name = param->name().size() + 3; 
  }

  for (unsigned i = 0; i < (unsigned)m_covarianceMatrix.GetNrows(); ++i ) {
    auto param = m_mps->at(i);
    if (param->name().find("_blind") != std::string::npos ) continue;
    double mean = param->mean(); 
    if ( param->isBlind() ) {
      double secretoffset =  m_mps->at(param->name()+"_blind")->mean();
      if (secretoffset == 0.) {
        INFO("\n\n\n Attempting to print a blind result!!! \n\n\n Skipping parameter for now, change this in FitResult.cpp");
        continue;
      }
      mean += secretoffset;
    }
    if( param->flag() == Flag::Free or
        param->flag() == Flag::Fix  or 
        param->flag() == Flag::Blind )  
    INFO( std::setw(longest_parameter_name)
          << param->name() << "     " << std::setw(5) << to_string<Flag>(param->flag())  
          << std::right << std::setw(13) << param->mean() << " ± "  
          << std::left  << std::setw(13) << (param->isFree() ?  m_mps->at(i)->err() : 0) 
          << " Pos err:" << std::setw(11) << (param->isFree() ? m_mps->at(i)->errPos() : 0) 
          << " Neg err:" << std::setw(11) << (param->isFree() ? m_mps->at(i)->errNeg() : 0)  
          << (param->isBlind() ? "  (BLIND)" : "") );
    }
  std::cout<<"\n"<<std::endl;
  */
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
    if ( ( param->isFree() || param->isBlind() || extended ) && param->err() > 1e-6 ) floating.push_back( param );
  }
  return floating;
}

TMatrixD FitResult::getReducedCovariance( const bool& extended ) const
{
  std::vector<unsigned int> floating_indices;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    auto param = m_mps->at(i);
    if ( ( param->isFree() || param->isBlind() || extended ) && param->err() > 1e-6 ) floating_indices.push_back( i );
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
double FitResult::Edm() const { return m_Edm; }
double FitResult::NCalls() const { return m_NCalls; }
int FitResult::status() const { return m_status; }
int FitResult::nParam() const { return m_nParam; }
int FitResult::nBins() const { return m_nBins; }

std::map<std::string, double> FitResult::observables() const { return m_observables; }
MinuitParameterSet* FitResult::mps() const { return m_mps; }

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

void FitResult::writeOptions( const std::string& output, const std::string& input)
{
  std::set<std::string> keys;
  std::ofstream output_stream( output );
  auto print_param = [&](auto& param) mutable {
    std::string rt= mysprintf( " %-7s %-12lf %-12lf", to_string<Flag>(param->flag()).c_str(), param->mean(), param->err() );
    if( param->minInit() != 0. && param->maxInit() != 0. ) rt += mysprintf(" %-12lf %-12lf", param->minInit(), param->maxInit()); 
    return rt; 
  };
  auto try_add_real    = [&](const auto& key) mutable
  {
    auto param = this->m_mps->find(key);
    if( param == nullptr ) return false; 
    output_stream << mysprintf("%-50s", param->name().c_str() ) << print_param(param) << "\n"; 
    keys.insert(key); 
    return true;  
  };
  auto try_add_complex = [&](const auto& key) mutable
  {
     auto param_re = this->m_mps->find(key+"_Re");
     auto param_im = this->m_mps->find(key+"_Im");
     if( param_re == nullptr || param_im == nullptr ) return false; 
     output_stream << mysprintf("%-50s", key.c_str()) << print_param(param_re) << print_param(param_im) << "\n"; 
     keys.insert(key+"_Re"); 
     keys.insert(key+"_Im"); 
     return true; 
  }; 

  if( input != "" )
  {
    processFile(input, [&](const std::string& line)
    {
      bool success = true; 
      auto tokens = split(line, ' ');
      if( tokens.size() == 0 ) { output_stream << line << "\n" ; return ;  } 
      auto key = tokens[0];
      if( ( m_mps->find(key) or m_mps->find(key +"_Re") or m_mps->find(key+"_Im") ) && tokens[1] != "=" )
      {
        success &= try_add_real(key) ? 1 : try_add_complex(key); 
      }
      else output_stream << line << "\n"; 
    }, '\0' );
  }
  else 
  {
    for( const auto& param : *m_mps )
    {
      auto name = param->name(); 
      if( keys.count(name) ) continue; 
      if( name.find("_Re") != std::string::npos && m_mps->find( replaceAll( name, "_Re", "_Im") ) ) 
        try_add_complex(replaceAll(name,"_Re","")); 
      else if( name.find("_Im") != std::string::npos && m_mps->find( replaceAll( name, "_Im", "_Re") ) ) 
        try_add_complex(replaceAll(name,"_Re",""));
      else try_add_real(name); 
    }
  }
}

