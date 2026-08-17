#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/MinuitExpression.h"
#include "AmpGen/ASTResolver.h"

#include <iomanip>

using namespace AmpGen;

namespace AmpGen {
  complete_enum(Flag, Free, Hide, Fix, CompileTimeConstant, Blind, Invalid)
}

MinuitParameter::MinuitParameter(const std::string &name, const Flag &fix, const double &mean, const double &step, const double &mi, const double &ma)
    : m_flag(fix), m_name(name), m_meanInit(mean), m_stepInit(step), m_minInit(mi), m_maxInit(ma) {
  DEBUG("Building parameter : " << name);
  resetToInit();
}

MinuitParameter::MinuitParameter(const std::string &name, const double &mean, const double &step, const double &mi, const double &ma)
    : MinuitParameter(name, Flag::Free, m_meanInit, m_stepInit, m_minInit, m_maxInit) {
  DEBUG("Building parameter : " << name);
  resetToInit();
}

MinuitParameter::~MinuitParameter() {
  for(auto &p : m_subscribers) p->m_parameter = nullptr;
}

Flag MinuitParameter::flag() const { return m_flag; }
double MinuitParameter::meanInit() const { return m_meanInit; }
double MinuitParameter::stepInit() const { return m_stepInit; }
double MinuitParameter::minInit() const { return m_minInit; }
double MinuitParameter::maxInit() const { return m_maxInit; }
double MinuitParameter::mean() const { return m_meanResult; }
double MinuitParameter::errPos() const { return m_errPosResult; }
double MinuitParameter::errNeg() const { return m_errNegResult; }
double MinuitParameter::err() const { return m_errResult; }
bool MinuitParameter::isFixed() const { return m_flag == Flag::Fix || m_flag == Flag::CompileTimeConstant; }
bool MinuitParameter::isFree() const { return m_flag == Flag::Free; }
bool MinuitParameter::isBlind() const { return m_flag == Flag::Blind; }

const std::string &MinuitParameter::name() const { return m_name; }

void MinuitParameter::fix() { m_flag = Flag::Fix; }
void MinuitParameter::scaleStep(const double &sf) {
  m_errResult *= sf;
  m_stepInit *= sf;
}
void MinuitParameter::setStepInit(const double &si) { m_stepInit = si; }

void MinuitParameter::setFree() {
  DEBUG("Setting parameter: " << m_name << " free");
  m_flag = Flag::Free;
}

void MinuitParameter::setVal(const double &val) {
  m_meanResult = val;
  double mu = mean();
  for(auto &subscriber : m_subscribers) subscriber->m_value = mu;
}

void MinuitParameter::setCurrentFitVal(double cfv) { setVal(cfv); }

void MinuitParameter::setInit(const double &val, const double &step) {
  m_meanInit = val;
  this->setVal(val);
  if(step != -1) m_stepInit = step;
}

void MinuitParameter::setResult(double fitMean, double fitErr, double fitErrNeg, double fitErrPos) {
  this->setVal(fitMean);
  m_errResult = fitErr;
  m_errPosResult = fitErrPos;
  m_errNegResult = fitErrNeg;
}

void MinuitParameter::setName(const std::string &name) { m_name = name; }

void MinuitParameter::broadcastToAll() {
  double nu = mean();
  for(auto &p : m_subscribers) p->m_value = nu;
}

void MinuitParameter::resetToInit() {
  this->setVal(m_meanInit);
  m_errResult = m_stepInit;
  m_errPosResult = -9999;
  m_errNegResult = -9999;
}

void MinuitParameter::setLimits(const double &min, const double &max) {
  m_minInit = min;
  m_maxInit = max;
}

std::ostream &AmpGen::operator<<(std::ostream &os, const MinuitParameter &par) {
  if(par.isBlind()) {
    return os << std::left << std::setw(60) << par.name() << " = " << std::right << std::setw(12) << " BLIND ± " << std::left << std::setw(12) << par.stepInit()
              << ((par.minInit() != 0 || par.maxInit() != 0) ? ("[" + std::to_string(par.minInit()) + ", " + std::to_string(par.maxInit())) + "]" : "")
              << " [flag=" << to_string<Flag>(par.flag()) << "]";
  } else {
    return os << std::left << std::setw(60) << par.name() << " = " << std::right << std::setw(12) << par.mean() << " ± " << std::left << std::setw(12)
              << par.stepInit()
              << ((par.minInit() != 0 || par.maxInit() != 0) ? ("[" + std::to_string(par.minInit()) + ", " + std::to_string(par.maxInit())) + "]" : "")
              << " [flag=" << to_string<Flag>(par.flag()) << "]";
  }
}

std::ostream &AmpGen::operator<<(std::ostream &os, const MinuitProxy &par) { return os << *par.parameter(); }

DEFINE_CAST(ExpressionParameter)

std::string ExpressionParameter::to_string(const ASTResolver *resolver) const {
  auto as_expression = dynamic_cast<const MinuitExpression *>(m_parameter.parameter());
  if(as_expression != nullptr) return as_expression->expression().to_string(resolver);

  if(resolver == nullptr and m_parameter.isValid()) return m_parameter->name();
  if(m_parameter.isValid() && m_parameter->flag() == Flag::CompileTimeConstant) return std::to_string(m_parameter->mean());
  return resolver->resolvedVariable(this);
}

std::string ExpressionParameter::name() const { return m_parameter->name(); }

void ExpressionParameter::resolve(ASTResolver &resolver) const {
  if(m_parameter.isValid()) {
    auto as_expression = dynamic_cast<const MinuitExpression *>(m_parameter.parameter());
    if(as_expression != nullptr) return as_expression->expression().resolve(resolver);
    if(m_parameter->flag() != Flag::CompileTimeConstant) resolver.resolve(*this);
  }
}

complex_t ExpressionParameter::operator()() const {
  if(!m_parameter.isValid()) {
    ERROR("Parameter does not have end-point");
    return complex_t(0., 0.);
  }
  return m_parameter->mean();
}
