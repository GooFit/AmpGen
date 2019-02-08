#ifndef AMPGEN_PHASESPACE_H
#define AMPGEN_PHASESPACE_H

#include <stddef.h>
#include <vector>

#include "AmpGen/Event.h"
#include "AmpGen/EventType.h"
#include "TRandom.h"

namespace AmpGen
{
  /**@class PhaseSpace
   * Phase-space generator taken from the ROOT routine, from Rene Brun and Valerio Filippini, 
   * which was originally based on the GENBOD routine of CERNLIB with
   * modifications such that unweighted events are generated,
   * which saves a large amount of CPU. Further modified to remove the dependencies on TLorentzVector 
   * and to only use the AmpGen::Event classes, and to add the option to include a time dependence.  
   */

  class PhaseSpace
  {
    public:
      PhaseSpace()  = default; 
      PhaseSpace( const EventType& type, TRandom* rand = nullptr );
      ~PhaseSpace() = default;
      
      bool setDecay( const double& m0, const std::vector<double>& mass );
      void setRandom( TRandom* rand ) { m_rand = rand; }
      size_t size() const { return m_nt; }
      AmpGen::Event makeEvent( const size_t& cacheSize = 0 );
      AmpGen::EventType eventType() const;

    private:
      unsigned int m_nt        = {0}; // number of decay particles
      double       m_mass[18]  = {0}; // masses of particles
      double       m_teCmTm    = {0}; // total energy in the C.M. minus the total mass
      double       m_wtMax     = {0}; // maximum weight 
      double       m_decayTime = {0}; // decay time
      TRandom*     m_rand      = {nullptr}; // Random number generator
      EventType    m_type;                  // EventType to generate

      double rndm() { return m_rand->Rndm(); }
      double q(double a, double b, double c) const;
  };
} // namespace AmpGen
#endif
