// @(#)root/physics:$Id$
// Author: Rene Brun , Valerio Filippini  06/09/2000

//_____________________________________________________________________________________
//
//  Utility class to generate n-body event,
//  with constant cross-section (default)
//  or with Fermi energy dependence (opt="Fermi").
//  The event is generated in the center-of-mass frame,
//  but the decay products are finally boosted
//  using the betas of the original particle.
//
//  The code is based on the GENBOD function (W515 from CERNLIB)
//  using the Raubold and Lynch method
//      F. James, Monte Carlo Phase Space, CERN 68-15 (1968)
//
// see example of use in $ROOTSYS/tutorials/physics/PhaseSpace.C
//
// Note that Momentum, Energy units are Gev/C, GeV

#include "AmpGen/PhaseSpace.h"

#include <RtypesCore.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "TRandom.h"
#include "AmpGen/Types.h"

const Int_t kMAXP = 18;

using namespace AmpGen;

double PhaseSpace::PDK( double a, double b, double c )
{
  double x = ( a - b - c ) * ( a + b + c ) * ( a - b + c ) * ( a + b - c );
  x        = sqrt( x ) / ( 2 * a );
  return x;
}

double PhaseSpace::Generate()
{
  std::array<double, kMAXP> rno;
  double pd[kMAXP];
  double invMas[kMAXP];

  rno[0] = 0;
  unsigned int n;

  double wt = m_wtMax;
  do {
    wt     = m_wtMax;
    rno[0] = 0;
    for ( n = 1; n < m_nt - 1; n++ ) rno[n] = rndm(); // m_nt-2 random numbers
    rno[m_nt - 1]                           = 1;
    std::sort( rno.begin() + 1, rno.begin() + m_nt );

    double sum = 0;
    for ( n = 0; n < m_nt; n++ ) {
      sum += m_mass[n];
      invMas[n] = rno[n] * m_teCmTm + sum;
    }
    for ( n = 0; n < m_nt - 1; n++ ) {
      pd[n] = PDK( invMas[n + 1], invMas[n], m_mass[n + 1] );
      wt *= pd[n];
    }
  } while ( wt < rndm() );

  m_decPro[0].SetPxPyPzE( 0, pd[0], 0, sqrt( pd[0] * pd[0] + m_mass[0] * m_mass[0] ) );

  unsigned int i = 1;
  unsigned int j;
  while ( 1 ) {
    m_decPro[i].SetPxPyPzE( 0, -pd[i - 1], 0, sqrt( pd[i - 1] * pd[i - 1] + m_mass[i] * m_mass[i] ) );
    double cZ   = 2 * rndm() - 1;
    double sZ   = sqrt( 1 - cZ * cZ );
    double angY = 2 * M_PI * rndm();
    double cY   = cos( angY );
    double sY   = sin( angY );
    for ( j = 0; j <= i; j++ ) {
      TLorentzVector& v = m_decPro[j];
      double x          = v.Px();
      double y          = v.Py();
      v.SetPx( cZ * x - sZ * y );
      v.SetPy( sZ * x + cZ * y ); // rotation around Z
      x        = v.Px();
      double z = v.Pz();
      v.SetPx( cY * x - sY * z );
      v.SetPz( sY * x + cY * z ); // rotation around Y
    }

    if ( i == ( m_nt - 1 ) ) break;
    double beta = pd[i] / sqrt( pd[i] * pd[i] + invMas[i] * invMas[i] );
    for ( j = 0; j <= i; j++ ) m_decPro[j].Boost( 0, beta, 0 );
    i++;
  }
  return wt;
}

TLorentzVector* PhaseSpace::GetDecay( const unsigned int& n )
{
  if ( n > m_nt ) return nullptr;
  return &( m_decPro[n] );
}

Bool_t PhaseSpace::SetDecay( TLorentzVector& P, const unsigned int& nt, const double* mass )
{
  std::vector<double> masses;
  masses.assign( mass, mass + nt );
  return SetDecay( P.Mag(), masses );
}

Bool_t PhaseSpace::SetDecay( const double& m0, const std::vector<double>& mass )
{
  unsigned int n;
  m_nt = mass.size();
  m_decPro.resize( m_nt, TLorentzVector( 0, 0, 0, 0 ) );

  m_teCmTm = m0; // total energy in C.M. minus the sum of the masses
  for ( n = 0; n < m_nt; n++ ) {
    m_mass[n] = mass[n];
    m_teCmTm -= mass[n];
  }

  if ( m_teCmTm <= 0 ) return 0 ; 

  double emmax = m_teCmTm + m_mass[0];
  double emmin = 0;
  double wtmax = 1;
  for ( n = 1; n < m_nt; n++ ) {
    emmin += m_mass[n - 1];
    emmax += m_mass[n];
    wtmax *= PDK( emmax, emmin, m_mass[n] );
  }
  m_wtMax = 1 / wtmax;

  return true;
}


PhaseSpace::PhaseSpace( const EventType& type, TRandom* rand ) : m_type( type ), m_rand( rand )
{
  SetDecay( type.motherMass(), type.masses() );
  if ( type.isTimeDependent() )
    m_decayTime = 6.582119514 / ( ParticlePropertiesList::get( type.mother() )->width() * pow( 10, 13 ) );
}

AmpGen::EventType PhaseSpace::eventType() const { return m_type; }

AmpGen::Event PhaseSpace::makeEvent( const unsigned int& cacheSize )
{
  Generate();
  AmpGen::Event newEvent( m_type.eventSize(), cacheSize );
  for ( unsigned int i = 0; i < m_nt; ++i ) {
    newEvent.set( i, {m_decPro[i].Px(), m_decPro[i].Py(), m_decPro[i].Pz(), m_decPro[i].E()} );
  }
  newEvent.setGenPdf( 1 );
  if ( m_type.isTimeDependent() ) newEvent.set( 4 * m_nt, m_rand->Exp( m_decayTime ) );
  return newEvent;
}
