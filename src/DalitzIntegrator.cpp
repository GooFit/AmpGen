#include "AmpGen/DalitzIntegrator.h"

#include <TLorentzVector.h>

#include "AmpGen/Kinematics.h"

using namespace AmpGen;

DalitzIntegrator::DalitzIntegrator( const double& mother, const double& m0, const double& m1, const double& m2 )
    : m_mother( mother ), m_massA( m0 ), m_massB( m1 ), m_massC( m2 ), m_rand( nullptr )
{
  setMin();
}

void DalitzIntegrator::set( const double& mother, const double& m0, const double& m1, const double& m2 )
{
  m_massA = m0;
  m_massB = m1;
  m_massC = m2;
  setMother( mother );
}
double DalitzIntegrator::sqDp1( const Event& evt ) const
{
  TLorentzVector p1( ( evt.address( 0 ) ) );
  TLorentzVector p2( ( evt.address( 4 ) ) );
  TLorentzVector p3( ( evt.address( 8 ) ) );
  TLorentzVector pA = p1 + p2;
  return acos( 2 * ( pA.Mag() - m_min ) / ( m_max - m_min ) - 1 ) / M_PI;
}
double DalitzIntegrator::sqDp2( const Event& evt ) const
{
  TLorentzVector p1( evt.address( 0 ) );
  TLorentzVector p2( evt.address( 4 ) );
  TLorentzVector p3( evt.address( 8 ) );
  TLorentzVector pA = p1 + p2;
  return acos( AmpGen::Product( p1, p3, pA ) / sqrt( AmpGen::Product( p1, p1, pA ) * AmpGen::Product( p3, p3, pA ) ) ) /
         M_PI;
}

void DalitzIntegrator::setEvent( const sqCo& x, double* event ) const
{
  double mAB = getMAB( x );
  double sA  = m_massA * m_massA;
  double sB  = m_massB * m_massB;
  double sC  = m_massC * m_massC;
  double pA  = sqrt( mAB * mAB / 4. - ( sA + sB ) / 2. + ( sA - sB ) * ( sA - sB ) / ( 4 * mAB * mAB ) );
  double eA  = sqrt( pA * pA + sA );
  double eB  = sqrt( pA * pA + sB );
  double eC  = ( m_mother * m_mother - sC - mAB * mAB ) / ( 2 * mAB );
  double pC  = sqrt( eC * eC - sC );

  event[2]  = pA;
  event[3]  = eA;
  event[6]  = -pA;
  event[7]  = eB;
  event[9]  = pC * sin( M_PI * x.second );
  event[10] = pC * cos( M_PI * x.second );
  event[11] = eC;
}

void DalitzIntegrator::setMin()
{
  m_min = m_massA + m_massB;
  m_max = m_mother - m_massC;
}

void DalitzIntegrator::setMother( const double& m )
{
  m_mother = m;
  setMin();
}

double DalitzIntegrator::getMAB( sqCo coords ) const
{
  double m = coords.first;
  return m_min + ( m_max - m_min ) * ( cos( M_PI * m ) + 1 ) / 2;
}

double DalitzIntegrator::J( const sqCo& coords ) const
{
  double mAB         = getMAB( coords );
  double sA          = m_massA * m_massA;
  double sB          = m_massB * m_massB;
  double sC          = m_massC * m_massC;
  double pA          = sqrt( mAB * mAB / 4. - ( sA + sB ) / 2. + ( sA - sB ) * ( sA - sB ) / ( 4 * mAB * mAB ) );
  double eC          = ( m_mother * m_mother - sC - mAB * mAB ) / ( 2 * mAB );
  double pC          = sqrt( eC * eC - sC );
  double mPrime_     = coords.first;
  double thetaPrime_ = coords.second;
  double deriv1      = ( M_PI / 2. ) * ( m_max - m_min ) * sin( M_PI * mPrime_ );
  double deriv2      = M_PI * sin( M_PI * thetaPrime_ );

  double j = 4 * deriv1 * deriv2 * pA * pC * mAB / ( 1000 * 1000 * 1000 );
  return j;
}
