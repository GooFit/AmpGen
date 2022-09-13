#include "AmpGen/PhaseSpace.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Types.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"

#include "TRandom.h"
const Int_t kMAXP = 18;

using namespace AmpGen;

PhaseSpace::PhaseSpace( const EventType& type, TRandom* rand ) : 
  m_rand(rand), m_type(type) 
{
  setDecay( type.motherMass(), type.masses() );
  if ( type.isTimeDependent() ){
    DEBUG("Generating with time-dependence");
    m_decayTime = ParticlePropertiesList::get( type.mother() )->lifetime();
  }
}

PhaseSpace::PhaseSpace( const Particle& particle, TRandom* rand ) : 
  PhaseSpace( particle.eventType(), rand ) {}

double PhaseSpace::q( double m, double m1, double m2 ) const
{
  return 0.5 * sqrt( m*m - 2*m1*m1 - 2*m2*m2 + (m1*m1-m2*m2)*(m1*m1-m2*m2)/(m*m) );
}

void PhaseSpace::fill( double* output, unsigned width  )
{
  std::array<double, kMAXP> rno;
  std::array<double, kMAXP> pd;
  std::array<double, kMAXP> invMas; 
  rno[0] = 0;
  size_t n;
  double wt = m_wtMax;
  do {
    wt     = m_wtMax;
    rno[0] = 0;
    for( n = 1; n < m_nt - 1; n++ ) rno[n] = m_rand->Rndm(); // m_nt-2 random numbers
    rno[m_nt - 1]                          = 1;
    std::sort( rno.begin() + 1, rno.begin() + m_nt );
    double sum = 0;
    for ( n = 0; n < m_nt; n++ ) {
      sum += m_mass[n];
      invMas[n] = rno[n] * m_teCmTm + sum;
    }
    for ( n = 0; n < m_nt - 1; n++ ) {
      pd[n] = q( invMas[n + 1], invMas[n], m_mass[n + 1] );
      wt *= pd[n];
    }
  } while ( wt < m_rand->Rndm() );
  output[0*width] = 0;
  output[1*width] = pd[0];
  output[2*width] = 0;
  output[3*width] = sqrt( pd[0] * pd[0] + m_mass[0] * m_mass[0] );
  for(size_t i = 1 ; i != m_nt ; ++i ){  
    output[(4*i+0)*width] = 0;
    output[(4*i+1)*width] = -pd[i-1];
    output[(4*i+2)*width] = 0;
    output[(4*i+3)*width] = sqrt( pd[i-1] * pd[i-1] + m_mass[i] * m_mass[i] );
    const auto cZ   = 2 * m_rand->Rndm() - 1;
    const auto sZ   = sqrt( 1 - cZ * cZ );
    const auto angY = 2 * M_PI * m_rand->Rndm();
    const auto cY   = cos(angY);
    const auto sY   = sin(angY);
    const auto beta  = (i == m_nt-1) ? 0 : pd[i] / sqrt( pd[i] * pd[i] + invMas[i] * invMas[i] );
    const auto gamma = (i == m_nt-1) ? 1 : 1./sqrt( 1 - beta*beta);    
    for (size_t j = 0; j <= i; j++ ) {
      const auto x          = output[(4*j+0)*width];
      const auto y          = output[(4*j+1)*width];
      const auto z          = output[(4*j+2)*width];
      const auto E          = output[(4*j+3)*width];
      output[(4*j+0)*width] = cY * (cZ * x - sZ * y ) - sY * z;
      output[(4*j+1)*width] = gamma*( sZ * x + cZ * y + beta * E );
      output[(4*j+2)*width] = sY * (cZ * x - sZ * y )  + cY * z;
      output[(4*j+3)*width] = gamma*( E + beta * (sZ *x + cZ*y) );
    }
  }
  if ( m_type.isTimeDependent() ) output[4 * m_nt * width] = m_rand->Exp( m_decayTime );
}

Event PhaseSpace::makeEvent()
{
  Event rt(4*m_nt + m_type.isTimeDependent());
  fill( rt.address()  );
  rt.setGenPdf( 1 ); 
  return rt;
}

bool PhaseSpace::setDecay( const double& m0, const std::vector<double>& mass )
{
  m_nt = mass.size();
  m_teCmTm = m0;
  for (size_t n = 0; n < m_nt; n++) {
    m_mass[n] = mass[n];
    m_teCmTm -= mass[n];
  }
  if ( m_teCmTm <= 0 ) return 0 ; 
  double emmax = m_teCmTm + m_mass[0];
  double emmin = 0;
  double wtmax = 1;
  for (size_t n = 1; n < m_nt; n++) {
    emmin += m_mass[n - 1];
    emmax += m_mass[n];
    wtmax *= q( emmax, emmin, m_mass[n] );
  }
  m_wtMax = 1 / wtmax;
  return true;
}

EventType PhaseSpace::eventType() const { return m_type; }

