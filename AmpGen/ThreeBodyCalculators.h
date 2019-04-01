#ifndef AMPGEN_THREEBODYCALCULATORS_H
#define AMPGEN_THREEBODYCALCULATORS_H

#include "AmpGen/DalitzIntegrator.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Tensor.h"

class TGraph;

namespace AmpGen
{
  class MinuitParameterSet;

  class ThreeBodyCalculator
  {
  private:
    struct PartialWidth {
      CoherentSum fcs;
      DalitzIntegrator integrator;
      CompiledExpression< std::complex<double>, const real_t*, const real_t* > totalWidth;
      EventType type;
      std::vector<CompiledExpression< std::complex<double>, const real_t*, const real_t*>> partialWidths; 
      double getWidth( const double& m );
      PartialWidth( const EventType& type, MinuitParameterSet& mps );
      Expression spinAverageMatrixElement( const std::vector<TransitionMatrix<std::complex<double>>>& elements,
                                           DebugSymbols* msym );
    };
    Expression calculateSAME( const std::string& particle );

    double       m_min;
    double       m_max;
    double       m_norm;
    double       m_step;
    size_t       m_nKnots;
    std::string  m_name;
    std::vector<PartialWidth> m_widths;
    MinuitParameterSet* m_mps; 
  public:
    ThreeBodyCalculator( const std::string& head, MinuitParameterSet& mps, const size_t& nKnots=999, const double& min=-1, const double& max=-1 );

    TGraph* widthGraph( const double& mNorm=-1 );
    TGraph* widthGraph( const size_t& steps, const double& min, const double& max );
    TGraph* runningMass( const double& mass, const double& min, const double& max, const size_t& nSteps, const size_t& nSubtractions=2 );
    TGraph* fastRunningMass( const double& mass, const double& min, const double& max, const size_t& nSteps, const size_t& nSubtractions=2 );

    double getWidth( const double& m );

    void updateRunningWidth( MinuitParameterSet& mps, const double& mNorm = 0 );    
    void setNorm( const double& mNorm );
    void setAxis( const size_t& nKnots, const double& min, const double& max );
    void prepare();
    void makePlots(const double& mass=-1, const size_t& x=0, const size_t& y=0);
    void debug( const double& m, const double& theta );
  };
} // namespace AmpGen
#endif /* end of include guard: AMPGEN_THREEBODYCALCULATORS_H */
