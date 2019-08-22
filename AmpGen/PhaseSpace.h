#ifndef AMPGEN_PHASESPACE_H
#define AMPGEN_PHASESPACE_H

#include <stddef.h>
#include <vector>

#include "AmpGen/Event.h"
#include "AmpGen/EventType.h"
#include "TRandom.h"

namespace AmpGen
{
  /** @class PhaseSpace
      Phase-space generator taken from the ROOT routine, from Rene Brun and Valerio Filippini, 
      which was originally based on the GENBOD routine of CERNLIB with
      modifications such that unweighted events are generated,
      which saves a large amount of CPU. Further modified to remove the dependencies on TLorentzVector 
      and to only use the AmpGen::Event classes, and to add the option to include a time dependence.  
   */
  class Particle;

  class PhaseSpace
  {
    public:
      PhaseSpace()  = default;                                    ///< Empty constructor 
      PhaseSpace(const EventType& type, TRandom* rand = gRandom); ///< Construct a phase space generator from an EventType 
      PhaseSpace(const Particle&  type, TRandom* rand = gRandom); ///< Construct a phase space generator from a Particle 
      
      bool setDecay( const double& m0, const std::vector<double>& mass ); ///< Set the parameters of this phase space generator
      void setRandom( TRandom* rand ) { m_rand = rand; } ///< Set the random number used by this phase space generator 
      size_t size() const { return m_nt; }               ///< Return the number of decay products
      Event makeEvent( const size_t& cacheSize = 0 );    ///< Make an event in this phase space. 
      EventType eventType() const;                       ///< Returns the EventType that this phase space is generating

    private:
      size_t       m_nt        = {0};                    ///< Number of particles in the final state
      double       m_mass[18]  = {0};                    ///< Masses of particles in the final state
      double       m_teCmTm    = {0};                    ///< Total energy in the rest frame minus the total mass
      double       m_wtMax     = {0};                    ///< Maximum weight of an event
      double       m_decayTime = {0};                    ///< Proper decay time of the particle
      TRandom*     m_rand      = {nullptr};              ///< Random number generator
      EventType    m_type;                               ///< EventType to generate

      double q(double a, double b, double c) const;      ///< Breakup momentum of a particle with mass a decaying to decay products with masses b and c.  
  };
} // namespace AmpGen
#endif
