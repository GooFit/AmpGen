#ifndef AMPGEN_MINIMISER_H
#define AMPGEN_MINIMISER_H

// Minimiser class using ROOT::Math::Minimiser
// allows generic use of Minuit1, Minuit2, other algorithms

#include <functional>
#include <iostream>
#include <vector>

#include "TMatrixTSym.h"

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

  class IExtendLikelihood;
  class MinuitParameter;
  class MinuitParameterSet;

  class Minimiser
  {
  private:
  template <typename T>
  struct HasGetVal
  {
    typedef char YesType[1];
    typedef char NoType[2]; 
    template <typename C> static YesType& test( decltype(&C::getVal) ) ;
    template <typename C> static NoType& test(...);
    enum { value = sizeof(test<T>(0)) == sizeof(YesType) };
  };
  
  public:
    template <typename TYPE> typename std::enable_if_t<HasGetVal<TYPE>::value, void> setFunction( TYPE& fcn )
    {
      m_theFunction = [&fcn]() { return fcn.getVal(); };
    }
    template <typename TYPE> typename std::enable_if_t<!HasGetVal<TYPE>::value, void> setFunction(TYPE& fcn)
    {
      m_theFunction = [&fcn](){ return fcn() ; } ;
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
    void addExtendedTerm( IExtendLikelihood* term );
    TMatrixTSym<double> covMatrix() const;
    TMatrixTSym<double> covMatrixFull() const;
    double operator()( const double* par );
    double FCN() const;
    MinuitParameterSet* parSet() const;
    int status() const;
    ROOT::Minuit2::Minuit2Minimizer* minimiserInternal();
    void GradientTest();

  
  private:
    MinuitParameterSet*         m_parSet       = {nullptr};
    std::function<double(void)> m_theFunction;
    ROOT::Minuit2::Minuit2Minimizer*  m_minimiser    = {nullptr};
    std::vector<double>         m_covMatrix    = {0};
    std::vector<unsigned int>   m_mapping      = {};
    int          m_status     = {0};
    unsigned int m_nParams    = {0};
    unsigned int m_printLevel = {0};
    double       m_ll_zero    = {0};
    bool         m_normalise  = {false};
    bool         m_debug      = {false};
    std::vector<IExtendLikelihood*> m_extendedTerms;
  };
} // namespace AmpGen
#endif
//
