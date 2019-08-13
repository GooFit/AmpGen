#ifndef AMPGEN_MINUITPARAMETER_H
#define AMPGEN_MINUITPARAMETER_H
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT

#include <iostream>
#include <string>
#include "AmpGen/enum.h"

namespace AmpGen
{
  class MinuitParameterSet;
  declare_enum( Flag, Free, Hide, Fix, CompileTimeConstant )

  class MinuitParameter
  {
  public: 

    MinuitParameter() = default;
    MinuitParameter(const std::string& name, const Flag& flag, const double& mean, const double& step,
                     const double& min = 0, const double& max = 0 );
    MinuitParameter(const std::string& name, const double& mean, const double& step,
                     const double& min = 0, const double& max = 0 );

    Flag flag() const;
    bool isFixed() const;
    bool isFree()  const;
    const std::string& name() const;

    double meanInit() const;
    double stepInit() const;
    double minInit() const;
    double maxInit() const;
    double err() const;
    double errPos() const;
    double errNeg() const;
    double* vp() { return &m_meanResult ; } 
    
    void setInit( const double& init );
    void setStepInit( const double& si );
    void setFree() ;
    void scaleStep( const double& sf );
    void fix();
    void setCurrentFitVal( double cfv );
    void setLimits( const double& min, const double& max );
    void setResult( double fitMean, double fitErr, double fitErrPos, double fitErrNeg );
    void resetToInit();
    void setName( const std::string& name );
    virtual double mean() const;
    virtual operator double() const { return m_meanResult; }
    virtual ~MinuitParameter() = default;
    
    friend class MinuitParameterSet;
  private:
    Flag m_flag;
    std::string m_name;
    double m_meanInit;
    double m_stepInit;
    double m_minInit;
    double m_maxInit;
    double m_meanResult;
    double m_errPosResult;
    double m_errNegResult;
    double m_errResult;
  };

  class MinuitProxy
  {
  public:
    void update() { m_value = m_parameter->mean(); }
    MinuitParameter* ptr() { return m_parameter; }
    operator double() const { return m_parameter->mean(); }
    double getFast() const { return m_value ; }
    MinuitProxy( MinuitParameter* param = nullptr ) : m_parameter( param ) { if( m_parameter != nullptr ) update(); }
    MinuitParameter* operator->() { return m_parameter; }
    const MinuitParameter* operator->() const { return m_parameter; }
  private:
    double m_value;
    MinuitParameter* m_parameter;
  };
  std::ostream& operator<<( std::ostream& os, const MinuitParameter& type );
} // namespace AmpGen

#endif
//
