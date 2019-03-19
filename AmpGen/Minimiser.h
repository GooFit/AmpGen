#ifndef AMPGEN_MINIMISER_H
#define AMPGEN_MINIMISER_H

// Minimiser class using ROOT::Math::Minimiser
// allows generic use of Minuit1, Minuit2, other algorithms

#include <functional>
#include <iostream>
#include <vector>

#include "TMatrixTSym.h"

namespace ROOT
{
  namespace Math
  {
    class Minimizer;
  }
}
class TGraph;

namespace AmpGen
{

  class IExtendLikelihood;
  class MinuitParameter;
  class MinuitParameterSet;

  class Minimiser
  {
  public:
    template <class TYPE>
    Minimiser( TYPE& fitFunction, MinuitParameterSet* mps ) : 
      m_parSet( mps )
    {
      m_theFunction = [&fitFunction]() { return fitFunction.getVal(); };
      prepare();
    }
    ~Minimiser() = default;
    
    unsigned int nPars() const;
    void prepare();
    bool doFit();
    void GradientTest();

    TGraph* scan( MinuitParameter* param, const double& min, const double& max, const double& step );

    void addExtendedTerm( IExtendLikelihood* term );
    TMatrixTSym<double> covMatrix() const;
    TMatrixTSym<double> covMatrixFull() const;
    double operator()( const double* par );
    double FCN() const;
    MinuitParameterSet* parSet() const;
    int status() const;
    ROOT::Math::Minimizer* minimiserInternal();
  
  private:
    void print( const double& LL );
    MinuitParameterSet* m_parSet;
    std::function<double( void )> m_theFunction;
    ROOT::Math::Minimizer* m_minimiser = {nullptr};
    std::vector<double> m_covMatrix    = {0};
    std::vector<unsigned int> m_mapping;
    std::vector<IExtendLikelihood*> m_extendedTerms;

    int m_status              = {0};
    unsigned int m_nParams    = {0};
    unsigned int m_lastPrint  = {0};
    unsigned int m_printLevel = {0};

  };
} // namespace AmpGen
#endif
//
