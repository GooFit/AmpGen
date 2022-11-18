#ifndef AMPGEN_TRANSITIONMATRIXELEMENT_H
#define AMPGEN_TRANSITIONMATRIXELEMENT_H

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Particle.h"

namespace AmpGen {
  class MatrixElement : public CompiledExpression<void(complex_v*, const size_t*, const real_t*, const real_v*)>
  {
    public:
      using amp_type    = CompiledExpression<void(complex_v*, const size_t*, const real_t*, const real_v*)>;
      MatrixElement() = default;
      MatrixElement(const Particle& dt, const TotalCoupling& coupling, const amp_type& amp); 
      MatrixElement(const Particle& dt, 
          const TotalCoupling& coupling, 
          const MinuitParameterSet& mps,
          const std::map<std::string, unsigned>& evtFormat, 
          const bool& debugThis=false);

      const std::vector<complex_v> operator()(const Event& event) const; 
      template <class... arg_types> auto operator()(arg_types... args ) const { return amp_type::operator()(args...) ; }
      void debug( const Event& event )    const; 
      const std::string decayDescriptor() const { return decayTree.decayDescriptor() ; }  
      Particle                                            decayTree;
      TotalCoupling                                       coupling;
      complex_t                                           coefficient;
      DebugSymbols                                        db; 
      bool                                                workToDo    = {false};
      unsigned                                            size        = {1}; // size of the returned object 
  };
  std::vector<size_t> processIndex(const std::vector<MatrixElement>& tm, const std::string& label);  
  size_t findIndex(const std::vector<MatrixElement>& tm, const std::string& decayDescriptor);
  std::vector<size_t> findIndices(const std::vector<MatrixElement>& tm, const std::string& decayDescriptor);
}
#endif
