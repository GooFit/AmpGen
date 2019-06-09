#include "AmpGen/DalitzIntegrator.h"

#include <TLorentzVector.h>
#include <Math/AdaptiveIntegratorMultiDim.h>
#include <Math/WrappedMultiTF1.h>
#include <iostream>

#include "AmpGen/Kinematics.h"
#include "AmpGen/Units.h"
#include "AmpGen/Event.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Projection.h"
#include "RtypesCore.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TRandom3.h"

using namespace AmpGen;
DalitzIntegrator::DalitzIntegrator( const double& s0, const double& s1, const double& s2, const double& s3 )
  : m_s0( s0 ), m_s1( s1 ), m_s2( s2 ), m_s3( s3 )
{
  DEBUG( "Setting up ISQDP with " << m_s0 << " " << m_s1 << "  "<< m_s2 << " "<< m_s3 );
  setMother(s0);
}

void DalitzIntegrator::set( const double& s0, const double& s1, const double& s2, const double& s3 )
{
  m_s1 = s1;
  m_s2 = s2;
  m_s3 = s3;
  setMother( s0 );
}
double DalitzIntegrator::sqDp1( const Event& evt ) const
{
  TLorentzVector p1( ( evt.address( 0 ) ) );
  TLorentzVector p2( ( evt.address( 4 ) ) );
  TLorentzVector p3( ( evt.address( 8 ) ) );
  TLorentzVector pA = p1 + p2;
  auto arg = 2 * ( pA.Mag() - m_min ) / ( m_max - m_min ) - 1;
  if( arg > 1 || arg < -1 ){
    ERROR("Argument: " << arg << " is out-of-bounds");
    return -1;
  }
  return acos(arg)/M_PI;
}
double DalitzIntegrator::sqDp2( const Event& evt ) const
{
  TLorentzVector p1( evt.address( 0 ) );
  TLorentzVector p2( evt.address( 4 ) );
  TLorentzVector p3( evt.address( 8 ) );
  TLorentzVector pA = p1 + p2;
  return acos( dotProduct( p1, p3, pA ) / sqrt( dotProduct( p1, p1, pA ) * dotProduct( p3, p3, pA ) ) ) / M_PI;
}

void DalitzIntegrator::setEvent( const sqCo& x, double* event ) const
{
  double mAB = getMAB( x );
  double pA  = safe_sqrt(0.25*mAB*mAB - 0.5*(m_s1 + m_s2) + (m_s1-m_s2)*(m_s1-m_s2)/(4*mAB*mAB));
  double eA  = sqrt(pA*pA + m_s1);
  double eB  = sqrt(pA*pA + m_s2);
  double eC  = ( m_s0 - m_s3 - mAB * mAB ) / ( 2 * mAB );
  double pC  = safe_sqrt( eC * eC - m_s3 );
  event[2]  = pA;
  event[3]  = eA;
  event[6]  = -pA;
  event[7]  = eB;
  event[9]  = pC * sin(M_PI * x.second);
  event[10] = pC * cos(M_PI * x.second);
  event[11] = eC;
}



void DalitzIntegrator::setMother( const double& s )
{
  m_s0 = s;
  m_min = sqrt(m_s1) + sqrt(m_s2); 
  m_max = sqrt(m_s0) - sqrt(m_s3);
}

double DalitzIntegrator::getMAB( sqCo coords ) const
{
  double m = coords.first;
  return m_min + ( m_max - m_min ) * ( cos( M_PI * m ) + 1 ) / 2;
}

double DalitzIntegrator::J( const sqCo& coords ) const
{
  double mAB         = getMAB( coords );
  double pA          = safe_sqrt( mAB * mAB / 4. - ( m_s1 + m_s2 ) / 2. + ( m_s1 - m_s2 ) * ( m_s1 - m_s2 ) / ( 4 * mAB * mAB ) );
  double eC          = ( m_s0 - m_s3 - mAB * mAB ) / ( 2 * mAB );
  double pC          = safe_sqrt( eC * eC - m_s3 );
  double mPrime_     = coords.first;
  double thetaPrime_ = coords.second;
  double deriv1      = ( M_PI / 2. ) * ( m_max - m_min ) * sin( M_PI * mPrime_ );
  double deriv2      = M_PI * sin( M_PI * thetaPrime_ );
  double j           = 4 * deriv1 * deriv2 * pA * pC * mAB / (GeV*GeV*GeV); 
  return j;
}

double DalitzIntegrator::J( const sqCo& coords, const double& s ) const
{
  double mAB         = getMAB( coords, s);
  double max         = sqrt(s) - sqrt(m_s3);
  double min         = sqrt(m_s1) + sqrt(m_s2);
  double pA          = sqrt( mAB*mAB/4. - (m_s1+m_s2)/2. + (m_s1-m_s2)*(m_s1-m_s2)/(4*mAB*mAB) );
  double eC          = ( m_s0 - m_s3 - mAB * mAB ) / (2*mAB);
  double pC          = safe_sqrt( eC * eC - m_s3 );
  double mPrime_     = coords.first;
  double thetaPrime_ = coords.second;
  double deriv1      = ( M_PI / 2. ) * ( max - min ) * sin( M_PI * mPrime_ );
  double deriv2      = M_PI * sin( M_PI * thetaPrime_ );
  double j           = 4 * deriv1 * deriv2 * pA * pC * mAB / (GeV*GeV*GeV); 
  return j;
}

double DalitzIntegrator::getMAB( sqCo coords, const double& s ) const
{
  double m = coords.first;
  double max         = sqrt(s) - sqrt(m_s3);
  double min         = sqrt(m_s1) + sqrt(m_s2);
  return min + ( max - min ) * ( cos( M_PI * m ) + 1 ) / 2;
}

void DalitzIntegrator::setEvent(const sqCo& x, double* event, const double& s) const
{
  double mAB = getMAB( x, s );
  double pA  = sqrt( mAB*mAB/4. - (m_s1+m_s2)/2. + (m_s1-m_s2)*(m_s1-m_s2)/(4*mAB*mAB) );
  double eA  = sqrt( pA * pA + m_s1 );
  double eB  = sqrt( pA * pA + m_s2 );
  double eC  = ( m_s0 - m_s3 - mAB * mAB ) / ( 2 * mAB );
  double pC  = sqrt( eC * eC - m_s3 );
  event[2]   = pA;
  event[3]   = eA;
  event[6]   = -pA;
  event[7]   = eB;
  event[9]   = pC * sin( M_PI * x.second );
  event[10]  = pC * cos( M_PI * x.second );
  event[11]  = eC;
}

void DalitzIntegrator::debug() const
{
  sqCo pos = { gRandom->Uniform(), gRandom->Uniform()};
  double event[12];
  for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
  setEvent( pos, event );
  INFO( "s0 = " << m_s0 << " " << sqrt(m_s0) <<  ", x = " << pos.first << " " << pos.second << " r = " << m_min << "  " << m_max );
  double mAB = getMAB(pos);
  double pA  = safe_sqrt( mAB * mAB / 4. - ( m_s1 + m_s2 ) / 2. + ( m_s1 - m_s2 ) * ( m_s1 - m_s2 ) / ( 4 * mAB * mAB ) );
  std::cout << "pA = " << pA << " mAB = " << mAB << std::endl; 
  for ( unsigned int i = 0; i < 12; ++i ) INFO( "Evt[" << i << "] = " << event[i] );
}

TH1D* DalitzIntegrator::makePlot( const std::function<double(const double*)>& fcn, const Projection& projection,
    const std::string& name, const size_t& nSamples )
{
  auto plot = projection.plot();
  double event[12];
  for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
  Event evtCache( 12 );
  for ( unsigned int i = 0; i < nSamples; ++i ) {
    sqCo pos = {gRandom->Uniform(), gRandom->Uniform()};
    setEvent( pos, event );
    evtCache.set( event );
    plot->Fill( projection( evtCache ), J( pos ) * fcn( event ) );
  }
  return plot;
}

TH2D* DalitzIntegrator::makePlot( const std::function<double(const double*)>& fcn, const Projection2D& projection,
    const std::string& name, const size_t& nSamples )
{
  auto plot = projection.plot();
  Event event( 12 );
  for ( unsigned int i = 0; i < nSamples; ++i ) {
    sqCo pos = {gRandom->Uniform(), gRandom->Uniform()};
    setEvent( pos, event );
    auto obs_cos = projection( event );
    plot->Fill( obs_cos.first, obs_cos.second, J( pos ) * fcn( event ) );
  }
  return plot;
}


double DalitzIntegrator::integrate_internal( TF2& fcn ) const
{
  ROOT::Math::WrappedMultiTF1 wf1( fcn );
  ROOT::Math::AdaptiveIntegratorMultiDim ig;
  ig.SetFunction( wf1 );
  ig.SetRelTolerance( 0.000001 );
  double xmin[] = {0, 0};
  double xmax[] = {1, 1};
  return ig.Integral(xmin,xmax);
}
