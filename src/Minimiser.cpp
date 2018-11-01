#include "AmpGen/Minimiser.h"

#include <TMatrixTBase.h>
#include <TMatrixTSym.h>
#include <iostream>
#include <string>

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TGraph.h"

using namespace AmpGen;

unsigned int Minimiser::nPars() const { return m_nParams; }

void Minimiser::print( const double& LL )
{
  INFO( "Iteration # " << m_lastPrint << "  FCN = " << round( LL, 3 ) << "  Edm = " << round( m_minimizer->Edm(), 3 )
                       << "  NCalls = " << m_minimizer->NCalls() );
}

double Minimiser::operator()( const double* xx )
{
  for ( unsigned int i = 0; i < m_mapping.size(); ++i ) {
    m_parSet->getParPtr( m_mapping[i] )->setCurrentFitVal( xx[i] );
  }
  double LL = m_theFunction();
  for ( auto& extendTerm : m_extendedTerms ) {
    LL -= 2 * extendTerm->getVal();
  }
  return LL;
}

void Minimiser::GradientTest()
{
  for ( unsigned int i = 0; i < m_mapping.size(); ++i ) {
    auto parameter = m_parSet->getParPtr( m_mapping[i] );
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
double Minimiser::FCN() const { return m_theFunction(); }

void Minimiser::prepare()
{
  std::string minimiser = NamedParameter<std::string>( "Minimiser::Minimiser", "Minuit2" );
  std::string algorithm = NamedParameter<std::string>( "Minimiser::Algorithm", "Hesse" );
  unsigned int maxCalls = NamedParameter<unsigned int>( "Minimiser::MaxCalls", 100000 );
  double tolerance      = NamedParameter<double>( "Minimiser::Tolerance", 0.01 );
  m_printLevel          = NamedParameter<unsigned int>( "Minimiser::PrintLevel", 4 );
  if ( m_minimizer != nullptr ) delete m_minimizer;
  m_minimizer = ROOT::Math::Factory::CreateMinimizer( minimiser, algorithm );
  DEBUG( "Error definition = " << m_minimizer->ErrorDef() );
  m_minimizer->SetMaxFunctionCalls( maxCalls );
  m_minimizer->SetMaxIterations( 100000 );
  m_minimizer->SetTolerance( tolerance );
  m_minimizer->SetPrintLevel( m_printLevel );

  m_mapping.clear();
  m_covMatrix.clear();

  unsigned int counter = 0;

  m_nParams = 0;

  DEBUG( "# parameters = " << m_parSet->size() );
  for ( unsigned int i = 0; i < m_parSet->size(); ++i ) {
    MinuitParameter* par = m_parSet->getParPtr( i );
    if ( par == nullptr ) {
      ERROR( "Parameter " << i << " doesn't exist!" );
      continue;
    }
    if ( par->iFixInit() != 0 ) continue;
    m_minimizer->SetVariable( counter, par->name(), par->mean(), par->stepInit() );
    if ( par->minInit() != 0 || par->maxInit() != 0 ) m_minimizer->SetVariableLimits( counter, par->minInit(), par->maxInit() );
    
    else if ( m_printLevel != 0 ) {
      INFO( "Parameter: " << std::left  << std::setw(60) << par->name() << " = " 
                          << std::right << std::setw(12) << par->mean() << " ± " 
                          << std::left  << std::setw(12) << par->stepInit() 
                          << ( ( par->minInit() != 0 || par->maxInit() != 0 ) ? 
                          ("["+ std::to_string(par->minInit()) + ", " + std::to_string(par->maxInit() ) ) : "" ) );
    }
    m_mapping.push_back( i );
    counter++;
  }
  m_nParams = counter;
  m_covMatrix.resize( m_nParams * m_nParams, 0 );
}

bool Minimiser::doFit()
{

  m_lastPrint = 999;
  ROOT::Math::Functor f( *this, m_nParams );
  for ( unsigned int i = 0; i < m_mapping.size(); ++i ) {
    MinuitParameter* par = m_parSet->getParPtr( m_mapping[i] );
    // INFO( "p["<<i<<"] = " << par->name() << "    " << par->mean() << "   " << par->err() );
    m_minimizer->SetVariable( i, par->name(), par->mean(), par->stepInit() );
    if ( par->minInit() != 0 || par->maxInit() != 0 ) {
      if ( m_printLevel != 0 )
        INFO( "Setting parameter limits : " << par->name() << " in [" << par->minInit() << ", " << par->maxInit()
                                            << "]" );
      m_minimizer->SetVariableLimits( i, par->minInit(), par->maxInit() );
    }
  }

  m_minimizer->SetFunction( f );
  m_minimizer->Minimize();
  for ( unsigned int i = 0; i < m_nParams; ++i ) {
    auto par = m_parSet->getParPtr( m_mapping[i] );
    ///    INFO( par->name() << " value=" << par->mean() << " gradient=" << (grad != 0 ? grad[i] : -1) );
    double error = *( m_minimizer->Errors() + i );
    par->setResult( *( m_minimizer->X() + i ), error, error, error );
    for ( unsigned int j = 0; j < m_nParams; ++j ) {
      m_covMatrix[i + m_nParams * j] = m_minimizer->CovMatrix( i, j );
    }
  }
  m_status = m_minimizer->Status();
  // if( m_printLevel != 0 )
  // INFO("Fit status = " << m_status );
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

void Minimiser::addExtendedTerm( IExtendLikelihood* m_term ) { m_extendedTerms.push_back( m_term ); }
