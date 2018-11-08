#ifndef AMPGEN_FITFRACTION_H
#define AMPGEN_FITFRACTION_H

#include <memory.h>
#include <memory>
#include <string>
#include <vector>

namespace AmpGen
{
  class EventType;
  class FastCoherentSum;
  class FastIncoherentSum;
  class Particle;

  struct FFCalculator {
    double operator()();
    std::vector<unsigned int> index_i;
    std::vector<unsigned int> index_j;
    std::vector<unsigned int> denom;
    FastCoherentSum* fcs;
    std::string name;
    FFCalculator( const std::string& _name, FastCoherentSum* _fcs, const std::vector<unsigned int> indices,
                  const std::vector<unsigned int> _denom );

    FFCalculator( const std::string& _name, FastCoherentSum* _fcs, const std::vector<unsigned int> _indexI,
                  const std::vector<unsigned int> _indexJ, const std::vector<unsigned int> _denom );
  };

  struct IFFCalculator {
    double operator()();
    unsigned int index;
    FastIncoherentSum* fcs;
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
