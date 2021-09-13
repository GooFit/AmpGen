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

#include "AmpGen/MetaUtils.h"
#include "AmpGen/enum.h"

/** @cond PRIVATE */
namespace ROOT
{
  namespace Minuit2
  {
    class Minuit2Minimizer;
  }
}
class TGraph;
/** @endcode */

namespace AmpGen
{
  declare_enum ( PrintLevel, Quiet, Info, Verbose, VeryVerbose ) 
  class ExtendLikelihoodBase;
  class MinuitParameter;
  class MinuitParameterSet;

  class Minimiser : public ROOT::Minuit2::MnTraceObject 
  {
  private:
    def_has_function(getVal)

  public:
    template <typename TYPE> void setFunction( TYPE& fcn )
    {
      if constexpr( has_getVal<TYPE>::value ) m_theFunction = [&fcn]() { return fcn.getVal(); };
      else m_theFunction = fcn;
    }

    template <typename TYPE> 
    Minimiser(TYPE& fitFunction, MinuitParameterSet* mps) : 
      m_parSet(mps)
    {
      setFunction(fitFunction);
      prepare();
    }
    
    Minimiser(std::function<double(void)>& fitFunction, MinuitParameterSet* mps) : 
      m_parSet(mps),
      m_theFunction(fitFunction)
    {
      prepare();
    }
    ~Minimiser() = default;
    
    unsigned int nPars() const;
    void prepare();
    void gradientTest();
    bool doFit();
    TGraph* scan( MinuitParameter* param, const double& min, const double& max, const double& step );
    void addExtendedTerm( ExtendLikelihoodBase* term );
    TMatrixTSym<double> covMatrix() const;
    TMatrixTSym<double> covMatrixFull() const;
    double operator()( const double* par );
    void operator()(int i, const ROOT::Minuit2::MinimumState & state)  override;
    double FCN() const;
    double Edm() const;
    double NCalls() const;  
    MinuitParameterSet* parSet() const;
    int status() const;
    ROOT::Minuit2::Minuit2Minimizer* minimiserInternal();
    void setPrintLevel( const PrintLevel& printLevel);
    void minos( MinuitParameter* param );
    ROOT::Fit::FitResult fitResult() const; 
  private:
    MinuitParameterSet*         m_parSet       = {nullptr};
    std::function<double(void)> m_theFunction;
    ROOT::Minuit2::Minuit2Minimizer*  m_minimiser    = {nullptr};
    std::vector<double>         m_covMatrix    = {0};
    std::vector<unsigned>       m_mapping      = {};
    int        m_status     = {0};
    unsigned   m_nParams    = {0};
    PrintLevel m_printLevel = {PrintLevel::Info};
    double     m_ll_zero    = {0};
    bool       m_normalise  = {false};
    std::vector<ExtendLikelihoodBase*> m_extendedTerms;

  };
} // namespace AmpGen
#endif
//
