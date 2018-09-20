// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

using namespace AmpGen;

MinuitParameter::Flag MinuitParameter::iFixInit() const { return m_iFixInit; }

double MinuitParameter::meanInit() const { return m_meanInit; }
double MinuitParameter::stepInit() const { return m_stepInit; }
double MinuitParameter::minInit() const { return m_minInit; }
double MinuitParameter::maxInit() const { return m_maxInit; }

double MinuitParameter::mean() const { return m_meanResult; }
double MinuitParameter::errPos() const { return m_errPosResult; }
double MinuitParameter::errNeg() const { return m_errNegResult; }
double MinuitParameter::err() const { return m_errResult; }
const std::string& MinuitParameter::name() const { return m_name; }

void MinuitParameter::setFree() { INFO("Setting parameter: " << m_name << " free" ) ; m_iFixInit = Flag::Float; }
bool MinuitParameter::hidden() const { return m_iFixInit == 1; }

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
  // ... in a format that can be read back in
  os << std::left << std::setw(65) << name() << "\t" << iFixInit() << "\t" << mean() << "\t" << err() << "\t" << minInit() << "\t" << maxInit();
}

MinuitParameter::MinuitParameter( const std::string& name, const MinuitParameter::Flag& fix, const double& mean, const double& step,
                                  const double& mi, const double& ma )
    : m_iFixInit( fix ), m_name( name ), m_meanInit( mean ), m_stepInit( step ), m_minInit( mi ), m_maxInit( ma )
{
  DEBUG( "Building parameter : " << name );

  resetToInit();
}

void MinuitParameter::resetToInit()
{
  m_meanResult = m_currentFitVal = m_meanInit;
  m_errResult                    = m_stepInit;
  m_errPosResult                 = -9999;
  m_errNegResult                 = -9999;
}
