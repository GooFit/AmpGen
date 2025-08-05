#include <ostream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "AmpGen/Spline.h"
#include "AmpGen/ASTResolver.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Property.h" 

using namespace AmpGen;

Spline::Spline(const std::string &name, const size_t &nKnots, const double &min, const double &max)
    : m_points(Parameter(name), 2 * nKnots), m_name(name), m_nKnots(nKnots), m_min(min), m_max(max) {}

Spline::Spline(const Spline &spline, const Expression &x, DebugSymbols *db)
    : m_points(spline.m_points), m_name(spline.m_name), m_nKnots(spline.m_nKnots), m_min(spline.m_min), m_max(spline.m_max), m_x(x), m_eval(eval(db)) {}

Expression Spline::operator()(const Expression &x, DebugSymbols *db) { return Spline(*this, x, db); }

SplineTransfer::~SplineTransfer() {
  if(m_transferMatrix != nullptr) gsl_matrix_free(m_transferMatrix);
}

Expression
AmpGen::getSpline(const std::string &name, const Expression &x, const std::string &arrayName, DebugSymbols *dbexpressions, const bool &continueSpline) {
  std::string spline_name = name + "::Spline::" + arrayName;
  using splineParam_t = std::tuple<unsigned, double, double>;
  splineParam_t splineParam = Property<splineParam_t>{nullptr, name + "::Spline", std::make_tuple(0, 0., 0.)};
  return Spline(spline_name, std::get<0>(splineParam), std::get<1>(splineParam), std::get<2>(splineParam))(x, dbexpressions);
}

Expression Spline::eval(DebugSymbols *db) const {
  Expression x = make_cse(m_x);
  double spacing = (m_max - m_min) / ((double)m_nKnots - 1.);
  Expression dx = Fmod(x - m_min, spacing);
  Expression bin = (x - m_min) / spacing;
  Expression continuedValue = m_points[m_nKnots - 1];

  Expression returnValue = Ternary(
    x > m_min && x < m_max,
    m_points[bin] + ((m_points[bin + 1] - m_points[bin]) / spacing - (m_points[bin + 1 + m_nKnots] + 2 * m_points[bin + m_nKnots]) * spacing / 6.) * dx
      + m_points[bin + m_nKnots] * dx * dx / 2. + dx * dx * dx * (m_points[bin + 1 + m_nKnots] - m_points[bin + m_nKnots]) / (6. * spacing),
    continuedValue);
  ADD_DEBUG(x, db);
  ADD_DEBUG(dx, db);
  ADD_DEBUG(bin, db);
  ADD_DEBUG(returnValue, db);
  ADD_DEBUG(m_points[bin], db);
  return make_cse(returnValue);
}

void SplineTransfer::print() const { INFO("Source: " << m_parameters[0]->name()); }

SplineTransfer::SplineTransfer(const size_t &address, const std::string &name, const unsigned int &N, const double &min, const double &max)
    : CacheTransfer(address, name), m_parameters(N, nullptr), m_nKnots(N), m_min(min), m_max(max) {
  unsigned int size = N - 2;
  gsl_matrix *M = gsl_matrix_alloc(size, size);
  gsl_matrix_set_all(M, 0.);
  m_transferMatrix = gsl_matrix_alloc(size, size);
  for(unsigned int i = 0; i < size; ++i) {
    gsl_matrix_set(M, i, i, 4.);
    if(i != size - 1) {
      gsl_matrix_set(M, i, i + 1, 1.);
      gsl_matrix_set(M, i + 1, i, 1.);
    }
  }
  gsl_permutation *p = gsl_permutation_alloc(size);
  int s;
  gsl_linalg_LU_decomp(M, p, &s);
  gsl_linalg_LU_invert(M, p, m_transferMatrix);
  gsl_permutation_free(p);
  gsl_matrix_free(M);
}

bool SplineTransfer::isConfigured() {
  return std::all_of(m_parameters.begin(), m_parameters.end(), [](auto &p) { return p != nullptr; });
}

void SplineTransfer::set(const unsigned int &N, MinuitParameter *f) { m_parameters[N] = f; }
void SplineTransfer::set(const unsigned int &N, const double &value) { m_parameters[N] = new MinuitParameter("dumb", Flag::Fix, value, 0); }

void SplineTransfer::transfer(CompiledExpressionBase *destination) {
  unsigned size = m_parameters.size() - 2;
  std::vector<double> mvectors(m_parameters.size(), 0);
  double spacing = (m_max - m_min) / double(m_parameters.size() - 1);

  for(unsigned i = 0; i != size; ++i) {
    for(unsigned j = 0; j != size; ++j) {
      double lj = m_parameters[j + 2]->mean() - 2 * m_parameters[j + 1]->mean() + m_parameters[j]->mean();
      mvectors[i + 1] += 6 * gsl_matrix_get(m_transferMatrix, i, j) * lj / (spacing * spacing);
    }
  }

  // CHECK_BLINDING
  for(unsigned i = 0; i < m_parameters.size(); ++i) {
    DEBUG("Knot GSL[" << i << "] : value = (" << m_parameters[i]->mean() << ", " << m_address + i << ") "
                      << " curvature = (" << round(mvectors[i], 10) << " , " << m_address + m_nKnots + i << ")");
    destination->setExternal(m_parameters[i]->mean(), m_address + i);
    destination->setExternal(mvectors[i], m_address + m_nKnots + i);
  }
}

void Spline::resolve(ASTResolver &resolver) const {
  resolver.resolve(*this);
  m_x.resolve(resolver);
  m_eval.resolve(resolver);
}

std::string Spline::to_string(const ASTResolver *resolver) const { return m_eval.to_string(resolver); }
Spline::operator Expression() { return Expression(std::make_shared<Spline>(*this)); }
complex_t Spline::operator()() const { return 0; }
