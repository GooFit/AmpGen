#ifndef AMPGEN_KMATRIX_H
#define AMPGEN_KMATRIX_H

#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"

namespace AmpGen
{
  struct poleConfig {
    Expression s;
    std::vector<Expression> couplings;
    std::vector<Expression> bl_factors;
    poleConfig(const Expression& s, const std::vector<Expression>& c = {})
        : s(s), couplings(c), bl_factors(c.size(), 1){};

    void add(const Expression& coupling, const Expression& bl_factor = 1)
    {
      couplings.push_back( coupling );
      bl_factors.push_back( bl_factor );
    }
    const Expression coupling(const size_t& i) const { return couplings[i] * bl_factors[i]; }
    const Expression g(const size_t& i) const { return couplings[i]; }
  };

  Tensor constructKMatrix(const Expression& s, const size_t& nChannels,
                          const std::vector<poleConfig>& poleConfigs);
  
  Expression phsp_twoBody( const Expression& s, const double& m0, const double& m1 );
  Expression phsp_fourPi( const Expression& s );
  Expression phsp_FOCUS( const Expression& s, const double& m0, const double& m1 );
  Expression gFromGamma( const Expression& m, const Expression& gamma, const Expression& rho );

  std::vector<Parameter> paramVector( const std::string& name, const unsigned int& nParam );

  Tensor getPropagator(const Tensor& kMatrix, const std::vector<Expression>& phaseSpace);


} // namespace AmpGen

#endif /* end of include guard: AMPGEN_KMATRIX_H */
