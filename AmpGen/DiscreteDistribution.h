#ifndef AMPGEN_DISCRETEDISTRIBUTION_H 
#define AMPGEN_DISCRETEDISTRIBUTION_H 1

#include <vector>
#include <stack> 
#include <tuple> 
#include "TRandom3.h"

// Generator of discrete distributions, 
// based on https://github.com/DavidPal/discrete-distribution 
// But modified to work with TRandom 
namespace AmpGen { 
  
  class DiscreteDistribution {
    public:
      DiscreteDistribution() = default; 
      DiscreteDistribution(const std::vector<double>& weights); 

      unsigned operator()(TRandom3* generator) const {
        const double number = generator->Uniform(); 
        size_t index = floor(m_buckets.size() * number);
        const auto& bucket = m_buckets[index];
        return number < std::get<2>(bucket) ? std::get<0>(bucket) : std::get<1>(bucket);
      }

      unsigned min() const; 
      unsigned max() const; 

      const std::vector<double>& probabilities() const; 
    private:
      void create_buckets(); 
      std::vector<double> m_probabilities;
      std::vector<std::tuple<unsigned, unsigned, double>> m_buckets;
  };
}

#endif
