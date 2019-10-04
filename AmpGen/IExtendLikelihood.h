#ifndef AMPGEN_IEXTENDLIKELIHOOD_H
#define AMPGEN_IEXTENDLIKELIHOOD_H

#include <string>
#include <vector>

namespace AmpGen
{
  class MinuitParameterSet;
  class MinuitParameter;
  class CoherentSum;

  class IExtendLikelihood
  {
  public:
    virtual ~IExtendLikelihood() = default;
    virtual double getVal() const = 0;
    virtual void configure( const std::string& configString, const AmpGen::CoherentSum& pdf,
                            const MinuitParameterSet& mps ) = 0;
    virtual IExtendLikelihood* create()                     = 0;
  };

  class GaussianConstraint : public IExtendLikelihood
  {
  public:
    double getVal() const override;
    void configure( const std::string& configString, const AmpGen::CoherentSum& pdf,
                    const MinuitParameterSet& mps ) override;
    IExtendLikelihood* create() override { return new GaussianConstraint(); }
    static std::string _id;

  private:
    MinuitParameter* m_param;
    double m_mean;
    double m_sigma;
  };

  class PartialWidthConstraint : public IExtendLikelihood
  {
  public:
    double getVal() const override;
    void configure( const std::string& configString, const AmpGen::CoherentSum& pdf,
                    const AmpGen::MinuitParameterSet& mps ) override;
    IExtendLikelihood* create() override { return new PartialWidthConstraint(); }
    static std::string _id;

  private:
    const CoherentSum* m_pdf;
    double m_ratio;
    double m_weight;
    std::vector<unsigned int> m_denComponents;
    std::vector<unsigned int> m_numComponents;
  };

  class LASSO : public IExtendLikelihood
  {
  public:
    double getVal() const override;
    void configure( const std::string& configString, const AmpGen::CoherentSum& pdf,
                    const MinuitParameterSet& mps ) override;
    IExtendLikelihood* create() override { return new LASSO(); }
    static std::string _id;

  private:
    double m_lambda;
    const CoherentSum* m_pdf;
  };
} // namespace AmpGen

#endif
