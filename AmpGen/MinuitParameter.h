#ifndef AMPGEN_MINUITPARAMETER_H
#define AMPGEN_MINUITPARAMETER_H

#include <iostream>
#include <string>
#include <set>
#include "AmpGen/enum.h"
#include "AmpGen/Expression.h"

namespace AmpGen {
  class MinuitParameterSet;
  class MinuitProxy;
  declare_enum(Flag, Free, Hide, Fix, CompileTimeConstant, Blind);

  class MinuitParameter {
  public:
    MinuitParameter() = default;
    MinuitParameter(const std::string &name, const Flag &flag, const double &mean, const double &step, const double &min = 0, const double &max = 0);
    MinuitParameter(const std::string &name, const double &mean, const double &step, const double &min = 0, const double &max = 0);
    virtual MinuitParameter clone() const { return *this; }
    virtual ~MinuitParameter();
    Flag flag() const;
    bool isFixed() const;
    bool isFree() const;
    bool isBlind() const;
    const std::string &name() const;

    double meanInit() const;
    double stepInit() const;
    double minInit() const;
    double maxInit() const;
    double err() const;
    double errPos() const;
    double errNeg() const;
    double *vp() { return &m_meanResult; }

    void setInit(const double &init, const double &step = -1);
    void setStepInit(const double &si);
    void setFree();
    void scaleStep(const double &sf);
    void fix();
    void setCurrentFitVal(double cfv);
    void setLimits(const double &min, const double &max);
    void setVal(const double &val);
    void setResult(double fitMean, double fitErr, double fitErrNeg, double fitErrPos);
    void resetToInit();
    void setName(const std::string &name);
    void broadcastToAll();
    virtual double mean() const;
    virtual operator double() const { return m_meanResult; }
    void setFromMinuitState(const double *x) {
      if(m_minuitIndex != -1) m_meanResult = x[m_minuitIndex];
    }
    void setMinuitIndex(const int &index) { m_minuitIndex = index; }
    int index() const { return m_minuitIndex; }

    friend class MinuitParameterSet;
    friend class MinuitProxy;

    void subscribe(MinuitProxy *sub) { m_subscribers.insert(sub); }
    void unsubscribe(MinuitProxy *sub) { m_subscribers.erase(sub); }
    void print_proxies() const {
      std::cout << m_name << "proxies: ";
      for(auto &p : m_subscribers) std::cout << p << " ";
      std::cout << std::endl;
    }

  protected:
    Flag m_flag;
    std::string m_name = {""};
    double m_meanInit = {0};
    double m_stepInit = {0};
    double m_minInit = {0};
    double m_maxInit = {0};
    double m_meanResult = {0};
    double m_errPosResult = {0};
    double m_errNegResult = {0};
    double m_errResult = {0};
    int m_minuitIndex = {-1};
    std::set<MinuitProxy *> m_subscribers;
  };

  class MinuitProxy {
  public:
    void update() {
      if(m_parameter != nullptr) m_value = m_parameter->mean();
    }
    MinuitProxy(const MinuitProxy &other) {
      m_parameter = other.m_parameter;
      m_value = other.m_value;
      if(m_parameter != nullptr) m_parameter->subscribe(this);
    }
    ~MinuitProxy() {
      if(m_parameter != nullptr) m_parameter->unsubscribe(this);
    }
    MinuitProxy &operator=(MinuitProxy &other) {
      m_parameter = other.m_parameter;
      m_value = other.m_value;
      if(m_parameter != nullptr) m_parameter->subscribe(this);
      return *this;
    }
    MinuitProxy &operator=(const MinuitProxy &other) {
      m_parameter = other.m_parameter;
      m_value = other.m_value;
      if(m_parameter != nullptr) m_parameter->subscribe(this);
      return *this;
    }
    MinuitProxy(MinuitProxy &&other) {
      m_parameter = other.m_parameter;
      m_value = other.m_value;
      if(m_parameter != nullptr) m_parameter->subscribe(this);
    }
    MinuitParameter *ptr() { return m_parameter; }
    operator double() const {
      // if( m_value != m_parameter->mean() ) WARNING("Desync'd proxy: " << m_parameter->name() << " " << this );
      return m_parameter == nullptr ? m_value : m_parameter->mean();
    }
    MinuitProxy(MinuitParameter *param = nullptr, const double &value = 0) : m_parameter(param), m_value(value) { update(); }
    MinuitParameter *operator->() { return m_parameter; }
    const MinuitParameter *operator->() const { return m_parameter; }
    MinuitParameter *parameter() { return m_parameter; }
    const MinuitParameter *parameter() const { return m_parameter; }
    bool isValid() const { return m_parameter != nullptr; }
    double m_value;

    friend class MinuitParameter;

  private:
    MinuitParameter *m_parameter{nullptr};
  };

  class ExpressionParameter : public IExpression {
  public:
    ExpressionParameter(const MinuitProxy &proxy) : m_parameter(proxy) {}
    std::string to_string(const ASTResolver *resolver = nullptr) const override;
    void resolve(ASTResolver &resolver) const override;
    complex_t operator()() const override;
    operator Expression() const;
    std::string name() const;

  private:
    MinuitProxy m_parameter;
  };

  std::ostream &operator<<(std::ostream &os, const MinuitParameter &);
  std::ostream &operator<<(std::ostream &os, const MinuitProxy &);
} // namespace AmpGen

#endif
//
