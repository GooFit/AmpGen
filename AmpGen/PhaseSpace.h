#ifndef AMPGEN_PHASESPACE_H
#define AMPGEN_PHASESPACE_H

// @(#)root/physics:$Id$
// Author: Rene Brun , Valerio Filippini  06/09/2000

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Phase Space Generator, based on the GENBOD routine of CERNLIB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "TLorentzVector.h"
#include "TRandom.h"

namespace AmpGen
{
  /**@class PhaseSpace
   * Phase-space generator based on ROOT routine, with
   * modifications such that unweighted events are generated,
   * which saves a large amount of CPU
   */

  class PhaseSpace
  {

  private:
    unsigned int m_nt; // number of decay particles
    double m_mass[18]; // masses of particles
    double m_teCmTm;   // total energy in the C.M. minus the total mass
    double m_wtMax;    // maximum weigth
    std::vector<TLorentzVector> m_decPro;
    EventType m_type;
    TRandom* m_rand;
    double m_decayTime;
    double rndm() { return m_rand->Rndm(); }
    double PDK( double a, double b, double c );

  public:
    PhaseSpace() : m_nt( 0 ), m_teCmTm( 0. ), m_wtMax( 0. ), m_rand( nullptr ) {}
    PhaseSpace( const PhaseSpace& gen );
    PhaseSpace& operator=( const PhaseSpace& gen );
    PhaseSpace( const EventType& type, TRandom* rand = nullptr );
    bool SetDecay( const double& m0, const std::vector<double>& mass );
    bool SetDecay( TLorentzVector& P, const unsigned int& nt, const double* mass );
    double Generate();
    size_t size() const { return m_nt; }
    TLorentzVector* GetDecay( const unsigned int& n );

    double GetWtMax() const { return m_wtMax; }

    AmpGen::Event makeEvent( const unsigned int& cacheSize = 0 );
    void setRandom( TRandom* rand ) { m_rand = rand; }
    AmpGen::EventType eventType() const;
  };
} // namespace AmpGen
#endif
