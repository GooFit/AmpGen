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
      DiscreteDistribution(const std::vector<double>& weights){
        normalize_weights(weights);
        create_buckets();
      }

      unsigned operator()(TRandom3* generator) const {
        const double number = generator->Uniform(); 
        size_t index = floor(m_buckets.size() * number);
        const auto& bucket = m_buckets[index];
        return number < std::get<2>(bucket) ? std::get<0>(bucket) : std::get<1>(bucket);
      }

      unsigned min() const { return 0;} 
      unsigned max() const { return m_probabilities.empty() ? 0 : m_probabilities.size() - 1; }

      const std::vector<double>& probabilities() const { return m_probabilities; }

    private:
      void normalize_weights(const std::vector<double>& weights) {
        const double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        m_probabilities.reserve(weights.size());
        for (auto weight : weights) m_probabilities.push_back(weight / sum);
      }

      void create_buckets() {
        auto topAndPop = [](auto& stack){ const auto s = stack.top(); stack.pop(); return s; };
        const size_t N = m_probabilities.size();
        if (N <= 0) {
          m_buckets.emplace_back(0, 0, 0.0);
          return;
        }
        std::stack<std::pair<double, size_t>> small, large; 
        // Split probabilities into small and large
        for (auto probability : m_probabilities) {
          ( (probability < 1.0/N) ? small : large ) . emplace( probability, small.size() + large.size() );  
        }

        m_buckets.reserve(N);

        unsigned i = 0;
        while (!small.empty() && !large.empty()) {
          const auto s = topAndPop(small); 
          const auto l = topAndPop(large); 
          // Create a mixed bucket
          m_buckets.emplace_back(s.second, l.second, s.first + static_cast<double>(i) / N);
          // Calculate the length of the left-over segment
          const double left_over = s.first + l.first - static_cast<double>(1) / N;
          // Re-insert the left-over segment
          ( left_over < (1.0 / N) ? small : large ) .emplace(left_over, l.second); 
          ++i;
        }

        // Create pure buckets
        while (!large.empty()) {
          const auto l = topAndPop(large);
          // The last argument is irrelevant as long it's not a NaN.
          m_buckets.emplace_back(l.second, l.second, 0.0);
        }

        // This loop can be executed only due to numerical inaccuracies.
        // TODO: Find an example when it actually happens.
        while (!small.empty()) {
          const auto s = topAndPop(small); 
          // The last argument is irrelevant as long it's not a NaN.
          m_buckets.emplace_back(s.second, s.second, 0.0);
        }
      }

      // List of probabilities
      std::vector<double> m_probabilities;
      std::vector<std::tuple<unsigned, unsigned, double>> m_buckets;
  };
}
