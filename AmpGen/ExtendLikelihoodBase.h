#ifndef AMPGEN_IEXTENDLIKELIHOOD_H
#define AMPGEN_IEXTENDLIKELIHOOD_H

#include <string>
#include <vector>

namespace AmpGen
{
  class MinuitParameterSet;
  class MinuitParameter;
  class CoherentSum;

  class ExtendLikelihoodBase
  {
  public:
    virtual ~ExtendLikelihoodBase() = default;
    virtual double operator()() const = 0;
    virtual void configure( const std::string& configString, const MinuitParameterSet& mps ) = 0;
    virtual ExtendLikelihoodBase* create()                     = 0;
  };

  class GaussianConstraint : public ExtendLikelihoodBase
  {
  public:
    double operator()() const override;
    GaussianConstraint() = default; 
    void configure( const std::string& configString, const MinuitParameterSet& mps ) override;
    ExtendLikelihoodBase* create() override { return new GaussianConstraint(); }
    static std::string _id;

  private:
    MinuitParameter* m_param;
    double m_mean;
    double m_sigma;
  };

  class PartialWidthConstraint : public ExtendLikelihoodBase
  {
  public:
    double operator()() const override;
    PartialWidthConstraint( const CoherentSum* pdf=nullptr) : m_pdf(pdf) {} 
    void configure( const std::string& configString, const AmpGen::MinuitParameterSet& mps ) override;
    ExtendLikelihoodBase* create() override { return new PartialWidthConstraint(); }
    static std::string _id;

  private:
    const CoherentSum* m_pdf;
    double m_ratio;
    double m_weight;
    std::vector<unsigned int> m_denComponents;
    std::vector<unsigned int> m_numComponents;
  };

  class LASSO : public ExtendLikelihoodBase
  {
  public:
    double operator()() const override;
    LASSO(const CoherentSum* pdf=nullptr) :m_pdf(pdf){};
    void configure( const std::string& configString, const MinuitParameterSet& mps ) override;
    ExtendLikelihoodBase* create() override { return new LASSO(); }
    static std::string _id;

  private:
    double m_lambda;
    const CoherentSum* m_pdf;
  };
} // namespace AmpGen

#endif
