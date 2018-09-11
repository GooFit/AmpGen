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
} // namespace ROOT
class TGraph;

namespace AmpGen
{

  class IExtendLikelihood;
  class MinuitParameter;
  class MinuitParameterSet;

  class Minimiser
  {
  private:
    void print( const double& LL );

  protected:
    MinuitParameterSet* m_parSet;
    std::function<double( void )> m_theFunction;
    ROOT::Math::Minimizer* m_minimizer;
    std::vector<double> m_covMatrix;
    std::vector<unsigned int> m_mapping;
    std::vector<IExtendLikelihood*> m_extendedTerms;

    int m_status;
    unsigned int m_nParams;
    unsigned int m_lastPrint;
    unsigned int m_printLevel;

  public:
    template <class TYPE>
    Minimiser( TYPE& fitFunction, MinuitParameterSet* mps )
        : m_parSet( mps ), m_minimizer( nullptr ), m_covMatrix( 0 ), m_status( 0 ), m_lastPrint( 0 )
    {
      m_theFunction = [&fitFunction]() { return fitFunction.getVal(); };
      prepare();
    }
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
    ROOT::Math::Minimizer* minimiserInternal() { return m_minimizer; }
  };
} // namespace AmpGen
#endif
//
