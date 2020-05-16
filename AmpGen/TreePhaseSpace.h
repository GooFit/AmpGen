#ifndef AMPGEN_TREEPHASESPACE_H
#define AMPGEN_TREEPHASESPACE_H

#include <memory.h>
#include <stddef.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include <random>

#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Particle.h"

#include <TRandom3.h>
#include <TLorentzVector.h>

namespace AmpGen
{
  class Particle;
  class Event; 
  /** @class TreePhaseSpace
    @brief Generator of events where the phase space is decomposed into a series of subtrees.  
    @decription Generates events using the decomposition of the phase space 
    into a series of two-body phase spaces, where the invariant mass of the two-body space 
    can be generated according to, for example, a simplified Breit-Wigner model. 
    This can then be used to generate a more complex model using accept-reject, which 
    allows for a much more efficient generation of narrow peaks. 

    Currently, each of the available channels is given the same weight relative to phase space, 
    ideally feedback could be given from the generator phase to focus on the more efficient channels, 
    i.e. those that have larger contributions to the full amplitude.   
    */
  class TreePhaseSpace
  {
    
    public:
      struct Vertex 
      {
        enum Type { BW, Flat, Stable, QuasiStable};
        Vertex() = default; 
        Vertex(const Particle& particle, const double& min);
        Vertex(const Particle& particle, const double& min, const double& max); 
        double p() const; 
        double weight() const; 
        double genPdf(const Event& event) const; 
        void generate();
        void print(const unsigned& offset = 0) const;
        void place(Event& event);
        Event event(const unsigned& eventSize);
        void generateFullEvent();
        void setRhoMax(); 
        void setRandom(TRandom3* rnd);
        static Vertex make(const Particle& particle, Vertex* parent = nullptr); 
        Particle    particle;
        double min      = {0};
        double max      = {0}; 
        double phiMin   = {0};
        double phiMax   = {0}; 
        Type   type     = {Type::BW};
        unsigned index  = {999};
        double bwMass   = {0};
        double bwWidth  = {0};
        double s        = {0};
        std::shared_ptr<Vertex> left    = {nullptr}; 
        std::shared_ptr<Vertex> right   = {nullptr};
        TRandom3* rand  = {nullptr};
        std::vector<unsigned> indices; 
        TLorentzVector mom;
        bool isMultiBody = {false}; 
        PhaseSpace phsp; /// multibody phase to resort to for non two-body decomposition;   
      };
      
      explicit TreePhaseSpace(const EventType& type);
      TreePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm = nullptr );
      TreePhaseSpace(const std::vector<Particle>& decayChains, const EventType& type, TRandom* rndm = nullptr);

      void setRandom( TRandom* rand );
      Event makeEvent();
      size_t size() const;
      EventType eventType() const ;
      double genPdf( const Event& event) const ; 
      const Vertex& operator[](const unsigned i) const { return m_top[i]; }      

      void provideEfficiencyReport(const std::vector<bool>& report);
    private:
      std::vector<Vertex>   m_top;
      TRandom3*             m_rand      = {nullptr}; 
      EventType             m_type;         ///< EventType to generate
      std::discrete_distribution<> m_dice;  ///< 
      std::vector<double>   m_weights; 
      std::vector<unsigned> m_generatorRecord;  
      std::mt19937          m_gen;
  };
} // namespace AmpGen

#endif
