#include "AmpGen/Minimiser.h"

#include <TMatrixTBase.h>
#include <TMatrixTSym.h>
#include <TGraph.h>
#include <iostream>
#include <string>
#include <iomanip>

#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MinimumParameters.h>
#include <Fit/FitConfig.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include "AmpGen/ExtendLikelihoodBase.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ProfileClock.h"

using namespace AmpGen;
using namespace ROOT;

namespace AmpGen {   
  complete_enum ( PrintLevel, Quiet, Info, Verbose, VeryVerbose); 
};

unsigned int Minimiser::nPars() const { return m_nParams; }

void Minimiser::operator()(int i, const ROOT::Minuit2::MinimumState & state) 
{
  if( m_printLevel == PrintLevel::Quiet ) return; 

  if( i >= 0)
  {
    INFO( "Iteration  #  "
        << std::setw(3) << i << " - FCN = " <<  std::setw(16) << state.Fval()
        << " Edm = " <<  std::setw(12) << state.Edm() << " NCalls = " << std::setw(6) << state.NFcn() );
    if (state.HasCovariance() )
      INFO("Error matrix change = " << state.Error().Dcovar() );
  }
  else {
    INFO( "Iteration  #  "
        << std::setprecision(13)
        << " - FCN = " <<  std::setw(16) << state.Fval()
        << " Edm = " <<  std::setw(12) << state.Edm() << " NCalls = " << std::setw(6) << state.NFcn() );
    if (state.HasCovariance() && m_printLevel >= PrintLevel::Verbose )
      INFO("Error matrix change = " << state.Error().Dcovar() );
  }
  auto gradient = state.Gradient().Vec().Data(); 
  // Print parameter values during the minimisation (blind params are shifted by the blinding offset)
  if( m_printLevel >= PrintLevel::Verbose )
  { 
    INFO( bold_on << "ID  Parameter                                           Value         Step          Gradient" << bold_off ); 
    for(size_t iii = 0 ; iii < m_mapping.size(); ++iii) {
      auto par = m_parSet->at(m_mapping[iii]);
      if ( ! (par->isFree() || par->isBlind() )  ) continue;      
      if (par->isBlind()){
        double secretoffset = m_parSet->at(par->name()+"_blind")->mean();
        if (!secretoffset){
          INFO("Attempting to print a blind parameter without secretoffset. Check Minimiser.cpp");
          continue;
        }	
        INFO( par->name()<< " " << par->mean() + secretoffset << "(BLIND)" ); 
      }
      else{
        double v = state.Vec().Data()[iii]; 
        if( par->minInit() != 0 || par->maxInit() != 0 ) v = par->minInit() + 0.5 * ( par->maxInit() - par->minInit() ) * ( std::sin( v ) + 1 );  
        INFO( mysprintf("%-3d %-50s % 2.4e   % 2.4e   % 2.4e", iii, par->name().c_str(), v, state.Gradient().Gstep()[iii], state.Gradient().Vec()[iii] ) ) ;
      }
    }
    std::cout << std::endl; 
  }
}


double Minimiser::operator()( const double* xx )
{
  ProfileClock callTime;
  for(size_t i = 0; i < m_mapping.size(); ++i ) {
    m_parSet->at( m_mapping[i] )->setCurrentFitVal( xx[i] );
  }
  double LL = m_theFunction() ;
  for ( auto& extendTerm : m_extendedTerms ) LL -= 2 * (*extendTerm)();
  callTime.stop();
  return LL - m_ll_zero;
}

double Minimiser::FCN() const { return m_theFunction(); }
double Minimiser::Edm() const { return m_minimiser->Edm(); }
double Minimiser::NCalls() const {return m_minimiser->NCalls(); }

void Minimiser::gradientTest()
{
  for (size_t i = 0; i < m_mapping.size(); ++i) {
    auto parameter = m_parSet->at( m_mapping[i] );
    double m       = parameter->mean();
    parameter->setCurrentFitVal( parameter->meanInit() + parameter->stepInit() );
    double vp = FCN();
    parameter->setCurrentFitVal( parameter->meanInit() - parameter->stepInit() );
    double vm = FCN();
    if ( parameter->stepInit() == 0 ) {
      WARNING( "Step of free parameter: " << parameter->name() << " = 0" );
    } else {
      INFO( " dF/d{" << parameter->name() << "} = " << ( vp - vm ) / ( 2 * parameter->stepInit() ) );
    }
    parameter->setCurrentFitVal( m );
  }
}

void Minimiser::prepare()
{
  std::string algorithm = NamedParameter<std::string>( "Minimiser::Algorithm", "Hesse");
  size_t maxCalls       = NamedParameter<size_t>( "Minimiser::MaxCalls"  , 100000);
  double tolerance      = NamedParameter<double>( "Minimiser::Tolerance" , 0.01);
  m_printLevel          = NamedParameter<PrintLevel>( "Minimiser::PrintLevel", PrintLevel::Info); 
  unsigned printLevelMinuit2 = NamedParameter<unsigned>("Minimiser::Minuit2MinimizerPrintLevel", m_printLevel == PrintLevel::VeryVerbose ? 3 : 0 );
  if( m_printLevel == PrintLevel::Invalid )
  {
    FATAL("Requested print level is not valid");
  }
  m_normalise           = NamedParameter<bool>(   "Minimiser::Normalise",false);
  if ( m_minimiser != nullptr ) delete m_minimiser;
  m_minimiser = new Minuit2::Minuit2Minimizer(algorithm.c_str() );
  DEBUG( "Error definition = " << m_minimiser->ErrorDef() );
  m_minimiser->SetMaxFunctionCalls( maxCalls );
  m_minimiser->SetMaxIterations( 100000 );
  m_minimiser->SetTolerance( tolerance );
  m_minimiser->SetStrategy( 3 );
  //  m_minimiser->SetPrecision(std::numeric_limits<double>::epsilon());
  m_minimiser->SetPrintLevel( printLevelMinuit2 ); // turn off minuit printing 
  m_mapping.clear();
  m_covMatrix.clear();

  /*
  if (m_printLevel >= PrintLevel::VeryVerbose && 
      std::any_of( m_parSet->begin(), m_parSet->end(), [](const auto& p ){ return p->isBlind(); } ) )
  {
    FATAL("Minimiser::PrintLevel is !=0, this is incompatible with having blind parameters.");
  }
  */
  for(size_t i = 0 ; i < m_parSet->size(); ++i)
  {
    auto par = m_parSet->at(i);
    if ( ! (par->isFree() || par->isBlind())  ) continue;
    m_minimiser->SetVariable(m_mapping.size(), par->name(), par->mean(), par->stepInit());
    if ( par->minInit() != 0 || par->maxInit() != 0 )
      m_minimiser->SetVariableLimits( m_mapping.size(), par->minInit(), par->maxInit() );
    m_mapping.push_back(i);
    if ( m_printLevel == PrintLevel::VeryVerbose )  INFO( *par );
  }
  m_nParams = m_mapping.size();
  m_covMatrix.resize( m_nParams * m_nParams, 0 );
}

std::string covMatrixStatusString( ROOT::Minuit2::Minuit2Minimizer* mini )
{
  if( mini->CovMatrixStatus() == -1 ) return "not available (inversion failed or Hesse failed)";
  if( mini->CovMatrixStatus() ==  0 ) return "available but not positive defined status";
  if( mini->CovMatrixStatus() ==  1 ) return "only approximate";
  if( mini->CovMatrixStatus() ==  2 ) return "is full matrix but forced pos def";
  if( mini->CovMatrixStatus() ==  3 ) return "is full and accurate"; 
  return "unknown";
} 

std::string minuitStatusString( ROOT::Minuit2::Minuit2Minimizer* mini )
{
  if( mini->Status() == 0 ) return "converged"; 
  if( mini->Status() == 1 ) return "forced positive definite";
  if( mini->Status() == 2 ) return "Hesse is invalid";
  if( mini->Status() == 3 ) return "Edm is above maximum";
  if( mini->Status() == 4 ) return "Reached max calls";
  if( mini->Status() == 5 ) return "unknown";
  return "unknown";
} 


bool Minimiser::doFit()
{
  if( m_normalise ) m_ll_zero = m_theFunction();
  ROOT::Math::Functor f( *this, m_nParams );
  for (size_t i = 0; i < m_mapping.size(); ++i ) {
    MinuitParameter* par = m_parSet->at( m_mapping[i] );
    m_minimiser->SetVariable( i, par->name(), par->mean(), par->stepInit() );
    if ( par->minInit() != 0 || par->maxInit() != 0 )
      m_minimiser->SetVariableLimits( i, par->minInit(), par->maxInit() );
  }
  m_minimiser->SetFunction( f );
  dynamic_cast< Minuit2::Minuit2Minimizer* >(m_minimiser)->SetTraceObject( *this );
  m_minimiser->Minimize();
  for (size_t i = 0; i < m_nParams; ++i ) {
    auto par = m_parSet->at( m_mapping[i] );
    double error = *( m_minimiser->Errors() + i );
    par->setResult( *( m_minimiser->X() + i ), error, error, error );
    for ( unsigned int j = 0; j < m_nParams; ++j ) {
      m_covMatrix[i + m_nParams * j] = m_minimiser->CovMatrix( i, j );
    }
  }
  m_status = m_minimiser->Status();

  INFO("Minuit2Minimize: " << minuitStatusString(m_minimiser) << ", covariance matrix: " << covMatrixStatusString( m_minimiser ) ); 
  INFO("Status = " << m_status ); 
  INFO("FVAL   = " << FCN() );
  INFO("Edm    = " << Edm() );
  INFO("Nfcn   = " << m_minimiser->NCalls() );
  if( m_status != 0 && m_printLevel != PrintLevel::VeryVerbose )
  {
    WARNING("Fit has not converged, some clues from Minuit2 may not be printed due to low verbosity level.");
    WARNING("Suggest using Minimiser::PrintLevel VeryVerbose or higher");
  }
  bool runMinos = NamedParameter<bool>("Minimiser::RunMinos",false);
  if(runMinos)
  {
    for( unsigned i = 0 ; i != m_nParams; ++i ){
      double low  = 0;
      double high = 0;
      int status  = 0;
      m_minimiser->GetMinosError(i, low, high, status);
      auto param = m_parSet->at( m_mapping[i] );
      param->setResult( *param, param->err(), low, high  );
    }
  }

  if( m_printLevel >= PrintLevel::Info )
  { 
    unsigned longest_parameter_name = 10;  
    for(const auto& param : *m_parSet )
    {
      if( param->name().size() > longest_parameter_name ) longest_parameter_name = param->name().size() + 3; 
    }
    for( const auto& param : *m_parSet )
    { 
      double mean = param->mean(); 
      if ( param->isBlind() ) {
        double secretoffset =  m_parSet->at(param->name()+"_blind")->mean();
        if (secretoffset == 0.) {
          INFO("\n\n\n Attempting to print a blind result!!! \n\n\n Skipping parameter for now, change this in Minimiser.cpp");
          continue;
        }
        mean += secretoffset;
      }
      if( param->flag() == Flag::Free or
          param->flag() == Flag::Fix  or 
          param->flag() == Flag::Blind ){
        if( runMinos ) INFO( std::setw(longest_parameter_name)
            << param->name() << "     " << std::setw(5) << to_string<Flag>(param->flag())  
            << std::right << std::setw(13) << mean << " ± "  
            << std::left  << std::setw(13) << (param->isFree() ?  param->err() : 0) 
            << " Pos err:" << std::setw(11) << (param->isFree() ? param->errPos() : 0) 
            << " Neg err:" << std::setw(11) << (param->isFree() ? param->errNeg() : 0) );
        else INFO( std::setw(longest_parameter_name)
            << param->name() << "     " << std::setw(5) << to_string<Flag>(param->flag())  
            << std::right << std::setw(13) << mean << " ± "  
            << std::left  << std::setw(13) << (param->isFree() ?  param->err() : 0) );
      }
    }
  }
  return 1;
}

TGraph* Minimiser::scan( MinuitParameter* param, const double& min, const double& max, const double& step )
{

  double secretoffset = 0.;
  if ( param->isBlind() ){ 
    secretoffset =  m_parSet->at( param->name()+"_blind")->mean();
  }

  param->fix();
  TGraph* rt = new TGraph();
  for ( double sv = min; sv < max; sv += step ) {
    param->setCurrentFitVal( sv - secretoffset);
    doFit();
    rt->SetPoint( rt->GetN(), sv, FCN() );
  }
  return rt;
}

int Minimiser::status() const { return m_status; }

TMatrixTSym<double> Minimiser::covMatrix() const
{
  unsigned int internalPars = nPars();
  TMatrixTSym<double> matrix( internalPars );
  for ( unsigned int i = 0; i < internalPars; i++ ) {
    for ( unsigned int j = 0; j < internalPars; j++ ) {
      matrix( i, j ) = m_covMatrix[i + j * internalPars];
    }
  }
  return matrix;
}

TMatrixTSym<double> Minimiser::covMatrixFull() const
{
  unsigned int internalPars = m_parSet->size();
  TMatrixTSym<double> matrix( internalPars );
  for ( unsigned int i = 0; i < m_mapping.size(); ++i ) {
    for ( unsigned int j = 0; j < m_mapping.size(); j++ ) {
      matrix( m_mapping[i], m_mapping[j] ) = m_covMatrix[i + j * m_mapping.size()];
    }
  }
  return matrix;
}

MinuitParameterSet* Minimiser::parSet() const { return m_parSet; }

void Minimiser::addExtendedTerm( ExtendLikelihoodBase* m_term )
{
  m_extendedTerms.push_back( m_term );
}

ROOT::Minuit2::Minuit2Minimizer* Minimiser::minimiserInternal() { return m_minimiser; }

void Minimiser::minos( MinuitParameter* parameter )
{ 

  if( m_minimiser == nullptr ) ERROR("No minimiser");

  ROOT::Math::Functor f( *this, m_nParams );
  m_minimiser->SetFunction( f );
  dynamic_cast< Minuit2::Minuit2Minimizer* >(m_minimiser)->SetTraceObject( *this );
  unsigned int index = 0; 
  for( ; index != m_nParams; ++index )
  {
    if( parameter->name() == m_parSet->at( m_mapping[index] )->name() ) break;
  }
  double v0 = parameter->mean();
  if( index == m_nParams ) ERROR( parameter->name() << " not amongst free parameters for fit");
  double low{0}, high{0};
  int status =0;
  std::vector<double> init_values ( m_minimiser->X(), m_minimiser->X() + m_nParams); 
  bool paramHasLimits = parameter->minInit() !=0 or parameter->maxInit() !=0;
  if( paramHasLimits != 0 && ( v0 - parameter->err() < parameter->minInit() ) ) return;  
  if( paramHasLimits != 0 && ( v0 + parameter->err() > parameter->maxInit() ) ) return;  
  
  m_minimiser->GetMinosError(index, low, high, status);
  parameter->setResult( v0, parameter->err(), low, high );
  
  for( int i = 0 ; i != m_nParams; ++i ) m_parSet->at( m_mapping[i] )->setCurrentFitVal( init_values[i] );

  if (parameter->isBlind())
  {
    double secretoffset = m_parSet->at(parameter->name()+"_blind")->mean();
    INFO( parameter->name() << " (BLIND) " << 
        parameter->mean()+secretoffset << "  - " << 
        parameter->errNeg() << " + " << 
        parameter->errPos() );
  }
  else  
  {
    INFO( parameter->name()   << 
        parameter->mean()   << "  - " << 
        parameter->errNeg() << " + " << 
        parameter->errPos() );
  }
}

void Minimiser::setPrintLevel( const PrintLevel& printLevel)
{ 
  m_printLevel = printLevel; 
  if (m_printLevel == PrintLevel::VeryVerbose ){
    /*
    for (const auto& param : *m_parSet){
      if ( param->isBlind() ) FATAL("Minimiser::PrintLevel is == VeryVerbose, incompatible with having any Blind parameter");
    }
    */
    //m_minimiser->SetPrintLevel( 2 );
  }
  
}

namespace AmpGen { 
  namespace detail { 
    class FitResult : public ROOT::Fit::FitResult {
      public:
        FitResult( const unsigned nParams ) : ROOT::Fit::FitResult(){ 
          fCovMatrix.resize( nParams  *  (1 + nParams ) / 2 ); 
          fErrors.resize(nParams);
          fParams.resize(nParams); 
          fParNames.resize(nParams); }
        void set( ROOT::Minuit2::Minuit2Minimizer* mini, int status, double fcn)
        {
          fNdf       = NPar();
          unsigned it= 0;
          for (unsigned i = 0; i < fNdf; ++i)
            for (unsigned j = 0; j <= i; ++j)
              fCovMatrix[it++] = mini->CovMatrix(i,j);

          fCovStatus = mini->CovMatrixStatus();
          fEdm       = mini->Edm();
          fVal       = fcn;
          fStatus    = status; 
          fErrors.assign( mini->Errors(), mini->Errors() + NPar() ); 
          fParams.assign( mini->X(), mini->X() + NPar() ); 
          fNCalls    = mini->NCalls();
          fValid     = fStatus == 0; 
        } 
        void setName( const unsigned int& index, const std::string& name )
        {
          fParNames[index] = name; 
        }
        void setAsymmError(const unsigned int& index, const double& low, const double& high )
        {
          fMinosErrors[index] = std::make_pair(low,high);
        }
        void blind(const unsigned int& index, const double& offset )
        {
          fParams[index] += offset;  
        }
    };
  }
}

ROOT::Fit::FitResult Minimiser::fitResult() const 
{
  detail::FitResult fr(m_nParams); 
  if( m_minimiser == nullptr ) return fr;
  fr.set(m_minimiser, status(), FCN()); 
  for( int i = 0 ; i != m_mapping.size(); ++i )
  {
    auto p = m_parSet->at ( m_mapping[i] );
    fr.setName( i, p->name() );
    if( p->errPos() != p->errNeg() ) fr.setAsymmError(i, p->errNeg(), p->errPos() );
    if( p->isBlind() ) fr.blind(i, m_parSet->at( p->name() + "_blind" )->mean() ); 
  }
  return ROOT::Fit::FitResult(fr); 
}
