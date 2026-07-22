#ifndef AMPGEN_IEXTENDLIKELIHOOD_H
#define AMPGEN_IEXTENDLIKELIHOOD_H

#include <string>
#include <vector>

#include "AmpGen/MinuitParameter.h"

namespace AmpGen {
  class MinuitParameterSet;
  class CoherentSum;

  class ExtendLikelihoodBase {
  public:
    virtual ~ExtendLikelihoodBase() = default;
    virtual double operator()() const = 0;
    virtual void configure(const std::string &configString, const MinuitParameterSet &mps) = 0;
    virtual ExtendLikelihoodBase *create() = 0;
  };

  class GaussianConstraint : public ExtendLikelihoodBase {
  public:
    double operator()() const override;
    GaussianConstraint() = default;
    void configure(const std::string &configString, const MinuitParameterSet &mps) override;
    ExtendLikelihoodBase *create() override { return new GaussianConstraint(); }
    static std::string _id;

  private:
    MinuitProxy m_param;
    double m_mean;
    double m_sigma;
  };

  class LASSO : public ExtendLikelihoodBase {
  public:
    double operator()() const override;
    LASSO(const CoherentSum *pdf = nullptr) : m_pdf(pdf) {};
    void configure(const std::string &configString, const MinuitParameterSet &mps) override;
    ExtendLikelihoodBase *create() override { return new LASSO(); }
    static std::string _id;

  private:
    double m_lambda;
    const CoherentSum *m_pdf;
  };
} // namespace AmpGen

#endif
