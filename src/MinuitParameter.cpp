// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT
#include "AmpGen/MinuitParameter.h"

#include <iomanip>

#include "AmpGen/MsgService.h"

using namespace AmpGen;

MinuitParameter::MinuitParameter( const std::string& name, const MinuitParameter::Flag& fix, const double& mean, const double& step,
    const double& mi, const double& ma )
  : m_iFixInit( fix ), m_name( name ), m_meanInit( mean ), m_stepInit( step ), m_minInit( mi ), m_maxInit( ma )
{
  DEBUG( "Building parameter : " << name );
  resetToInit();
}

MinuitParameter::MinuitParameter( const std::string& name, const double& mean, const double& step, const double& mi, const double& ma )
  : MinuitParameter(name, MinuitParameter::Flag::Float, m_meanInit, m_stepInit, m_minInit, m_maxInit ) 
{
  DEBUG( "Building parameter : " << name );
  resetToInit();
}

MinuitParameter::Flag MinuitParameter::iFixInit() const { return m_iFixInit; }

double MinuitParameter::meanInit()                const { return m_meanInit; }
double MinuitParameter::stepInit()                const { return m_stepInit; }
double MinuitParameter::minInit()                 const { return m_minInit; }
double MinuitParameter::maxInit()                 const { return m_maxInit; }
double MinuitParameter::mean()                    const { return m_meanResult; }
double MinuitParameter::errPos()                  const { return m_errPosResult; }
double MinuitParameter::errNeg()                  const { return m_errNegResult; }
double MinuitParameter::err()                     const { return m_errResult; }

bool MinuitParameter::fixed() { return m_iFixInit == Flag::Fix; }

void MinuitParameter::fix() { m_iFixInit = Flag::Fix; }
void MinuitParameter::scaleStep( const double& sf )
{
  m_errResult *= sf;
  m_stepInit *= sf;
}
void MinuitParameter::setStepInit( const double& si )
{
  m_stepInit = si;
}

const std::string& MinuitParameter::name()        const { return m_name; }
bool MinuitParameter::hidden()                    const { return m_iFixInit == 1; }

void MinuitParameter::setFree() { INFO("Setting parameter: " << m_name << " free" ) ; m_iFixInit = Flag::Float; }

void MinuitParameter::setCurrentFitVal( double cfv ) { m_meanResult = cfv; }

void MinuitParameter::setInit( const double& val ) { m_meanInit = val; }

void MinuitParameter::setResult( double fitMean, double fitErr, double fitErrPos, double fitErrNeg )
{
  m_meanResult   = fitMean;
  m_errResult    = fitErr;
  m_errPosResult = fitErrPos;
  m_errNegResult = fitErrNeg;
}

void MinuitParameter::print( std::ostream& os ) const
{
  os << std::left << std::setw(65) << name() << "\t" << iFixInit() << "\t" << mean() << "\t" << err() << "\t" << minInit() << "\t" << maxInit();
}

void MinuitParameter::setName( const std::string& name ) { m_name = name; }

void MinuitParameter::resetToInit()
{ 
  m_meanResult   = m_meanInit;
  m_errResult    = m_stepInit;
  m_errPosResult = -9999;
  m_errNegResult = -9999;
}

void MinuitParameter::setLimits( const double& min, const double& max )
{
  m_minInit = min;
  m_maxInit = max;
}
