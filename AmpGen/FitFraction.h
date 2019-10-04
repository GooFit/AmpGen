#ifndef AMPGEN_FITFRACTION_H
#define AMPGEN_FITFRACTION_H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include "AmpGen/Types.h"
#include "AmpGen/ErrorPropagator.h"

namespace AmpGen
{
  class EventType;
  class Particle;

  class FitFraction
  {
    public:
      FitFraction(const std::string& line);
      FitFraction(const std::string& name, const double& frac, const double& err);
      FitFraction() = default;

      void setFracErr(const double& f, const double& e);
      double val() const;
      double err() const;
      std::string name() const;
      std::shared_ptr<Particle> particle() const;

    private:
      std::string m_name;
      double m_value;
      double m_error;
  };
  
  bool operator  <(const FitFraction& lhs, const FitFraction& rhs);
  bool operator  >(const FitFraction& lhs, const FitFraction& rhs);
  bool operator ==(const FitFraction& lhs, const FitFraction& rhs);
  std::ostream& operator<<(std::ostream& os, const FitFraction& obj);
  
  template <class pdf_type>
  struct FitFractionCalculator {
    struct fcalc {
      std::string name; 
      std::vector<size_t> i;
      std::vector<size_t> j;
      fcalc(const std::string& name, const std::vector<size_t>& i) :
        name(name), i(i), j(i) {};
      fcalc(const std::string& name, const std::vector<size_t>& i, const std::vector<size_t>& j) :
        name(name), i(i), j(j) {};
    };
    pdf_type* pdf;
    std::vector<fcalc> calculators;
    std::vector<size_t> normSet; 
    bool recalculateIntegrals;
    template <class... ARGS> void emplace_back( ARGS&&... args ){ calculators.emplace_back(args...); }
  
    FitFractionCalculator(pdf_type* pdf, 
                          const std::vector<size_t>& normSet, 
                          const bool& recalculateIntegrals=false) : 
      pdf(pdf),
      normSet(normSet),
      recalculateIntegrals(recalculateIntegrals) {}
    std::vector<double> operator()(){
      if ( recalculateIntegrals ) pdf->prepare();
      else pdf->transferParameters();
      std::vector<double> rv;
      double sum = 0; 
      for (size_t i = 0; i != calculators.size(); ++i){
        auto v = getVal(i);
        rv.push_back(v);
        sum += v;
      }
      rv.push_back(sum);
      return rv; 
    }
    real_t norm() const {
      complex_t sum = 0;
      for ( auto& i : normSet ) {
        for ( auto& j : normSet ) {
          sum += (*pdf)[i].coefficient * std::conj( (*pdf)[j].coefficient ) * ( j >= i ? pdf->norm(i, j) : std::conj(pdf->norm(j,i)) );
        }
      }
      return std::real(sum);
    }
    real_t getVal( const size_t& index, const bool& getImaginaryPart = false) const {
      complex_t sum = 0; 
      for ( auto& i : calculators[index].i ) {
        for ( auto& j : calculators[index].j ) {
          sum += (*pdf)[i].coefficient * std::conj( (*pdf)[j].coefficient ) * ( j >= i ? pdf->norm(i, j) : std::conj(pdf->norm(j,i)) );
        }
      }
      return (getImaginaryPart ? std::imag(sum) : std::real(sum) ) / norm(); 
    }
    std::vector<FitFraction> operator()(const std::string& name, const LinearErrorPropagator& linProp )
    {
      auto values = (*this)();
      auto errors = linProp.getVectorError(*this,calculators.size()+1);
      std::vector<FitFraction> fractions; 
      for(size_t i = 0; i < calculators.size(); ++i) 
        fractions.emplace_back(calculators[i].name, values[i],errors[i]);
      std::sort(fractions.begin(), fractions.end());
      std::reverse(fractions.begin(), fractions.end());
      fractions.emplace_back("Sum_"+name, *values.rbegin(), *errors.rbegin() );
      return fractions;
    }   
  }; 
} // namespace AmpGen

#endif
