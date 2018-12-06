#ifndef AMPGEN_FITFRACTION_H
#define AMPGEN_FITFRACTION_H

#include <memory.h>
#include <memory>
#include <string>
#include <vector>

namespace AmpGen
{
  class EventType;
  class CoherentSum;
  class IncoherentSum;
  class Particle;

  struct FFCalculator {
    double operator()();
    std::vector<size_t> index_i;
    std::vector<size_t> index_j;
    std::vector<size_t> denom;
    CoherentSum* fcs;
    std::string name;
    FFCalculator( const std::string& name, CoherentSum* fcs, const std::vector<size_t>& indices,
                  const std::vector<size_t>& denom );

    FFCalculator( const std::string& name, CoherentSum* fcs, const std::vector<size_t>& indexI,
                  const std::vector<size_t>& indexJ, const std::vector<size_t>& denom );
  };

  struct IFFCalculator {
    double operator()();
    size_t index;
    IncoherentSum* fcs;
    IFFCalculator( const size_t& index, IncoherentSum* fcs );
  };

  class FitFraction
  {
    std::string m_name;
    double m_value;
    double m_error;

  public:
    FitFraction( const std::string& line, const AmpGen::EventType& evtType );
    FitFraction( const std::string& name, const double& frac, const double& err );
    FitFraction() = default;

    void setFracErr( const double& f, const double& e )
    {
      m_value = f;
      m_error = e;
    }
    double val() const { return m_value; }
    double err() const { return m_error; }
    std::string name() const { return m_name; }
    std::shared_ptr<Particle> particle() const;
  };
} // namespace AmpGen

#endif
