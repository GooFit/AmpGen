#ifndef AMPGEN_MINUITPARAMETER_H
#define AMPGEN_MINUITPARAMETER_H
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT

#include <iostream>
#include <string>

namespace AmpGen
{

  class MinuitParameterSet;

  class MinuitParameter
  {
  public:
    enum Flag {
      Float, Hide, Fix, CompileTimeConstant
    };
  private:
    Flag m_iFixInit;
    std::string m_name;
    double m_meanInit;
    double m_stepInit;
    double m_minInit;
    double m_maxInit;
    double m_meanResult;
    double m_errPosResult;
    double m_errNegResult;
    double m_errResult;
    double m_currentFitVal;

  public:

    MinuitParameter() = default;
    MinuitParameter( const std::string& name, const Flag& iFixInit, const double& mean, const double& step,
                     const double& min = 0, const double& max = 0 );

    Flag iFixInit() const;
    bool hidden() const;
    const std::string& name() const;

    void scaleStep( const double& sf )
    {
      m_errResult *= sf;
      m_stepInit *= sf;
    }
    double meanInit() const;
    double stepInit() const;
    double minInit() const;
    double maxInit() const;

    void setInit( const double& init );
    void setStepInit( const double& si ){ m_stepInit = si; }
    void setFree() ;
    void fix() { m_iFixInit = Flag::Fix; }
    bool fixed() { return m_iFixInit == Flag::Fix; }
    virtual double mean() const;
    double err() const;
    double errPos() const;
    double errNeg() const;

    void setCurrentFitVal( double cfv );
    void setLimits( const double& min, const double& max )
    {
      m_minInit = min;
      m_maxInit = max;
    }
    void setResult( double fitMean, double fitErr, double fitErrPos, double fitErrNeg );
    void resetToInit();

    void print( std::ostream& os = std::cout ) const;

    virtual ~MinuitParameter() = default;
    void setName( const std::string& name ) { m_name = name; }
    friend class MinuitParameterSet;

    virtual operator double() const { return m_currentFitVal; }
  };

  class MinuitProxy
  {
    double m_value;
    MinuitParameter* m_parameter;

  public:
    void update() { m_value = m_parameter->mean(); }
    MinuitParameter* ptr() { return m_parameter; }
    //    operator double() const { return m_value ; }
    operator double() const { return m_parameter->mean(); }
    double getFast() const { return m_value ; }
    MinuitProxy( MinuitParameter* param = nullptr ) : m_parameter( param ) { update(); }
    MinuitParameter* operator->() { return m_parameter; }
    const MinuitParameter* operator->() const { return m_parameter; }
  };
} // namespace AmpGen

#endif
//
