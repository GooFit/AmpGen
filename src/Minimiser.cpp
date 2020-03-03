#include "AmpGen/Minimiser.h"

#include <TMatrixTBase.h>
#include <TMatrixTSym.h>
#include <TGraph.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <Minuit2/Minuit2Minimizer.h>

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

using namespace AmpGen;
using namespace ROOT;

unsigned int Minimiser::nPars() const { return m_nParams; }

double Minimiser::operator()( const double* xx )
{
  for(size_t i = 0; i < m_mapping.size(); ++i ) {
    m_parSet->at( m_mapping[i] )->setCurrentFitVal( xx[i] );
  }
  double LL = m_theFunction() ;
  for ( auto& extendTerm : m_extendedTerms ) {
    LL -= 2 * extendTerm->getVal();
  }
  return LL - m_ll_zero;
}

double Minimiser::FCN() const { return m_theFunction(); }

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
  double tolerance      = NamedParameter<double>( "Minimiser::Tolerance" , 1);
  m_printLevel          = NamedParameter<size_t>( "Minimiser::PrintLevel", 4);
  m_normalise           = NamedParameter<bool>(   "Minimiser::Normalise",false);
  if ( m_minimiser != nullptr ) delete m_minimiser;
  m_minimiser = new Minuit2::Minuit2Minimizer(algorithm.c_str() );
  DEBUG( "Error definition = " << m_minimiser->ErrorDef() );
  m_minimiser->SetMaxFunctionCalls( maxCalls );
  m_minimiser->SetMaxIterations( 100000 );
  m_minimiser->SetTolerance( tolerance );
  m_minimiser->SetPrintLevel( m_printLevel );
  m_mapping.clear();
  m_covMatrix.clear();
  for(size_t i = 0 ; i < m_parSet->size(); ++i)
  {
    auto par = m_parSet->at(i);
    if ( ! par->isFree() ) continue;
    m_minimiser->SetVariable(m_mapping.size(), par->name(), par->mean(), par->stepInit());
    if ( par->minInit() != 0 || par->maxInit() != 0 )
      m_minimiser->SetVariableLimits( m_mapping.size(), par->minInit(), par->maxInit() );
    m_mapping.push_back(i);
    if ( m_printLevel != 0 ) INFO( *par );
  }
  m_nParams = m_mapping.size();
  m_covMatrix.resize( m_nParams * m_nParams, 0 );
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
  if(NamedParameter<bool>("Minimiser::RunMinos",true)){
    for( unsigned i = 0 ; i != m_nParams; ++i ){
      double low  = 0;
      double high = 0;
      int status  = 0;
      m_minimiser->GetMinosError(i, low, high, status);
      auto param = m_parSet->at( m_mapping[i] );
      param->setResult( *param, param->err(), low, high  );
    }
  }

  if(NamedParameter<bool>("Minimiser::PrintResults",true))
    print_parameters_and_covariance();

  return 1;
}

TGraph* Minimiser::scan( MinuitParameter* param, const double& min, const double& max, const double& step )
{
  param->fix();
  TGraph* rt = new TGraph();
  for ( double sv = min; sv < max; sv += step ) {
    param->setCurrentFitVal( sv );
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

void Minimiser::addExtendedTerm( IExtendLikelihood* m_term )
{
  m_extendedTerms.push_back( m_term );
}

ROOT::Minuit2::Minuit2Minimizer* Minimiser::minimiserInternal() { return m_minimiser; }

void Minimiser::print_parameters_and_covariance() const {
  //print results from floating fitparameters
  INFO(TString::Format("%-2s %-60s  %-11s %-11s %-11s","#","parametername","value","error high","error low"));
  for( unsigned i = 0 ; i != m_nParams; ++i ){
    auto const param = m_parSet->at( m_mapping[i] );
    //printing the parameter-number allows to identify the parameter in the covariance matrix
    INFO(TString::Format("%-2u %-60s  %11.5g %11.5g %11.5g",i,param->name().data(),param->mean(),param->errPos(),param->errNeg()));
  }

  INFO(std::string(128, '*'));
  INFO(TString::Format("%-127s","* Correlation matrix: "));
  INFO(std::string(128, '*'));

  for (int i = -2; i < static_cast<int>(m_nParams); ++i ) {
    TString line = "";
    for(int j = -1; j < static_cast<int>(m_nParams); j++){
      if(i == -2 && j == -1) line += "par # |";
      else if (i == -2) line += TString::Format("   %-4d",j);
      else if(i == -1 && j == -1) line += "______|";
      else if(j == -1) line += TString::Format("%-6d|",i);
      else if(i == -1) line += "_______";
      else{
        auto const current_corr = m_minimiser->Correlation(m_mapping[i],m_mapping[j]);
        if(current_corr > 0.2 && current_corr < 0.5)line += TString::Format(" \033[1;33m%-+.3f\033[0m",current_corr);
        else if(current_corr > 0.5 && current_corr < 1.f)line += TString::Format(" \033[1;31m%-+.3f\033[0m",current_corr);
        else if(current_corr < -0.2 && current_corr > -0.5)line += TString::Format(" \033[1;34m%-+.3f\033[0m",current_corr);
        else if(current_corr < -0.5)line += TString::Format(" \033[1;35m%-+.3f\033[0m",current_corr);
        else line += TString::Format(" %-+.3f",current_corr);
      }
    }
    INFO(line);
  }

  INFO("One more time for copying into the options file");
  for(auto const& p_idx : m_mapping)
    INFO(TString::Format("%-60s  0 %-11.5g %-11.5g",m_parSet->at(p_idx)->name().data(),m_parSet->at(p_idx)->mean(),m_parSet->at(p_idx)->errPos()));

}
