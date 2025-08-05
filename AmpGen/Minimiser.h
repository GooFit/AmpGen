#ifndef AMPGEN_MINIMISER_H
#define AMPGEN_MINIMISER_H

// Minimiser class using ROOT::Math::Minimiser
// allows generic use of Minuit1, Minuit2, other algorithms

#include <functional>
#include <iostream>
#include <vector>

#include <TMatrixTSym.h>
#include <Fit/FitResult.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnTraceObject.h>
#include <Math/IFunction.h>

#include "AmpGen/MetaUtils.h"
#include "AmpGen/enum.h"
#include "AmpGen/Property.h"
#include "AmpGen/Configurable.h"

/** @cond PRIVATE */
namespace ROOT {
  namespace Minuit2 {
    class Minuit2Minimizer;
  }
}
class TGraph;
/** @endcode */

namespace AmpGen {
  make_enum(PrintLevel, Quiet, Info, Verbose, VeryVerbose);

  class ExtendLikelihoodBase;
  class MinuitParameter;
  class MinuitParameterSet;

  class Minimiser : public Configurable<Minimiser>, ROOT::Minuit2::MnTraceObject {
  private:
    def_has_function(getVal);
    def_has_function(grad);

      public : 
    
    Minimiser() = default; 
    ~Minimiser() = default;

    template <typename TYPE> void setFunction(TYPE &fcn) {
      if constexpr(has_getVal<TYPE>::value)
        m_theFunction = [&fcn]() { return fcn.getVal(); };
      else {
        m_theFunction = fcn;
      }
      if constexpr(std::is_convertible<TYPE *, ROOT::Math::IGradientFunctionMultiDimTempl<double> *>::value) m_fcnWithGrad = &fcn;
    }
    template <typename TYPE> Minimiser(TYPE &fitFunction, MinuitParameterSet *mps) : m_parSet(mps) {
      setFunction(fitFunction);
      prepare();
    }

    Minimiser(std::function<double(void)> &fitFunction, MinuitParameterSet *mps) : m_parSet(mps), m_theFunction(fitFunction) { prepare(); }

    unsigned int nPars() const;
    void prepare();
    void gradientTest();
    bool doFit();
    TGraph *scan(MinuitParameter *param, const double &min, const double &max, const double &step);
    void addExtendedTerm(ExtendLikelihoodBase *term);
    TMatrixTSym<double> covMatrix() const;
    TMatrixTSym<double> covMatrixFull() const;
    double operator()(const double *par);
    void operator()(int i, const ROOT::Minuit2::MinimumState &state) override;
    double FCN() const;
    double Edm() const;
    double NCalls() const;
    MinuitParameterSet *parSet() const;
    int status() const;
    void setPrintLevel(const PrintLevel &printLevel);
    void minos(MinuitParameter *param);
    ROOT::Fit::FitResult fitResult() const;

  private:
    using GradFcn = ROOT::Math::IGradientFunctionMultiDimTempl<double>;
    std::vector<ExtendLikelihoodBase *> m_extendedTerms;
    GradFcn *m_fcnWithGrad{nullptr};
    MinuitParameterSet *m_parSet{nullptr};
    std::function<double(void)> m_theFunction{nullptr};
    ROOT::Math::Minimizer *m_minimiser{nullptr};
    std::vector<double> m_covMatrix{0};
    std::vector<unsigned> m_mapping{};
    int m_status{0};
    unsigned m_nParams{0};

    Property<std::string> m_minimiserTool{this, "Minimiser::Minimiser", "Minuit2"};
    Property<std::string> m_algorithm{this, "Minimiser::Algorithm", "Migrad"};
    Property<unsigned> m_maxCalls{this, "Minimiser::MaxCalls", 100000};
    Property<double> m_tolerance{this, "Minimiser::Tolerance", 1.0};
    Property<PrintLevel> m_printLevel{this, "Minimiser::PrintLevel", PrintLevel::Info};
    Property<double> m_precision{this, "Minimiser::Precision", 1e-15};
    Property<unsigned> m_printLevelMinuit2{this, "Minimiser::Minuit2MinimizerPrintLevel", m_printLevel == PrintLevel::VeryVerbose ? 3u : 0u};
    Property<bool> m_runMinos{this, "Minimiser::RunMinos", false};
  };
} // namespace AmpGen
#endif
//
