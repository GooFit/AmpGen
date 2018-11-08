#ifndef AMPGEN_PHASESPACE_H
#define AMPGEN_PHASESPACE_H

// @(#)root/physics:$Id$
// Author: Rene Brun , Valerio Filippini  06/09/2000

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Phase Space Generator, based on the GENBOD routine of CERNLIB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stddef.h>
#include <vector>

#include "AmpGen/Event.h"
#include "AmpGen/EventType.h"
#include "TLorentzVector.h"
#include "TRandom.h"

namespace AmpGen
{
  /**@class PhaseSpace
   * Phase-space generator taken from the ROOT routine, from Rene Brun and Valerio Filippini, which was originally based on the GENBOD routine of CERNLIB with
   * modifications such that unweighted events are generated,
   * which saves a large amount of CPU
   */

  class PhaseSpace
  {

  private:
    unsigned int m_nt  = {0}; // number of decay particles
    double m_mass[18]  = {0}; // masses of particles
    double m_teCmTm    = {0}; // total energy in the C.M. minus the total mass
    double m_wtMax     = {0}; // maximum weight 
    double m_decayTime = {0}; // decay time
    EventType m_type;
    TRandom*  m_rand = {nullptr};
    
    std::vector<TLorentzVector> m_decPro;
    double rndm() { return m_rand->Rndm(); }
    double PDK( double a, double b, double c );

  public:
    PhaseSpace() = default; 
    PhaseSpace( const EventType& type, TRandom* rand = nullptr );
    bool SetDecay( const double& m0, const std::vector<double>& mass );
    bool SetDecay( TLorentzVector& P, const unsigned int& nt, const double* mass );
    double Generate();
    double GetWtMax() const { return m_wtMax; }
    void setRandom( TRandom* rand ) { m_rand = rand; }
    
    size_t size() const { return m_nt; }
    TLorentzVector* GetDecay( const unsigned int& n );
    AmpGen::Event makeEvent( const unsigned int& cacheSize = 0 );
    AmpGen::EventType eventType() const;
  };
} // namespace AmpGen
#endif
