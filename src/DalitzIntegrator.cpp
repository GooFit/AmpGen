#include "AmpGen/DalitzIntegrator.h"

#include <TLorentzVector.h>
#include <iostream>
#include <cmath>

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
#include "TGraph.h"
#include <math.h>

#if ENABLE_AVX
#define SQRT sqrt
#else 
#define SQRT std::sqrt
#endif

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
  TLorentzVector p1( evt.address( 0 ) );
  TLorentzVector p2( evt.address( 4 ) );
  TLorentzVector p3( evt.address( 8 ) );
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
  return std::acos( dotProduct( p1, p3, pA ) / std::sqrt( dotProduct( p1, p1, pA ) * dotProduct( p3, p3, pA ) ) ) / M_PI;
}

void DalitzIntegrator::setEvent( const sqCo& x, real_v* event ) const
{
  auto mAB = getMAB( x );
  auto pA  = safe_sqrt(0.25*mAB*mAB - 0.5*(m_s1 + m_s2) + (m_s1-m_s2)*(m_s1-m_s2)/(4.*mAB*mAB));
  auto eA  = SQRT(pA*pA + m_s1);
  auto eB  = SQRT(pA*pA + m_s2);
  auto eC  = ( m_s0 - m_s3 - mAB * mAB ) / ( 2. * mAB );
  auto pC  = safe_sqrt( eC * eC - m_s3 );
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
  using namespace std;
  m_min = std::sqrt(m_s1) + std::sqrt(m_s2); 
  m_max = std::sqrt(m_s0) - std::sqrt(m_s3);
}

real_v DalitzIntegrator::getMAB( const sqCo& coords ) const
{
  return m_min + ( m_max - m_min ) * ( cos( M_PI * coords.first ) + 1 ) / 2.;
}

real_v DalitzIntegrator::J( const sqCo& coords ) const
{
  auto mAB         = getMAB( coords );
  auto pA          = safe_sqrt( mAB * mAB / 4. - ( m_s1 + m_s2 ) / 2. + ( m_s1 - m_s2 ) * ( m_s1 - m_s2 ) / ( 4. * mAB * mAB ) );
  auto eC          = ( m_s0 - m_s3 - mAB * mAB ) / ( 2. * mAB );
  auto pC          = safe_sqrt( eC * eC - m_s3 );
  auto mPrime_     = coords.first;
  auto thetaPrime_ = coords.second;
  auto deriv1      = ( M_PI / 2. ) * ( m_max - m_min ) * sin( M_PI * mPrime_ );
  auto deriv2      = M_PI * sin( M_PI * thetaPrime_ );
  auto j           = 4. * deriv1 * deriv2 * pA * pC * mAB;
  return j;
}

real_v DalitzIntegrator::J( const sqCo& coords, const double& s ) const
{
  using namespace std; 
  auto mAB         = getMAB( coords, s);
  auto max         = std::sqrt(s)    - std::sqrt(m_s3);
  auto min         = std::sqrt(m_s1) + std::sqrt(m_s2);
  auto pA          = SQRT( mAB*mAB/4. - (m_s1+m_s2)/2. + (m_s1-m_s2)*(m_s1-m_s2)/(4*mAB*mAB) );
  auto eC          = ( m_s0 - m_s3 - mAB * mAB ) / (2*mAB);
  auto pC          = safe_sqrt( eC * eC - m_s3 );
  auto mPrime_     = coords.first;
  auto thetaPrime_ = coords.second;
  auto deriv1      = ( M_PI / 2. ) * ( max - min ) * sin( M_PI * mPrime_ );
  auto deriv2      = M_PI * sin( M_PI * thetaPrime_ );
  return 4 * deriv1 * deriv2 * pA * pC * mAB;
}

real_v DalitzIntegrator::getMAB( const sqCo& coords, const double& s ) const
{
  auto m           = coords.first;
  auto max         = std::sqrt(s) - std::sqrt(m_s3);
  auto min         = std::sqrt(m_s1) + std::sqrt(m_s2);
  return min + ( max - min ) * ( cos( M_PI * m ) + 1 ) / 2;
}

void DalitzIntegrator::setEvent(const sqCo& x, real_v* event, const double& s) const
{
  auto mAB = getMAB( x, s );
  auto pA  = SQRT( mAB*mAB/4. - (m_s1+m_s2)/2. + (m_s1-m_s2)*(m_s1-m_s2)/(4*mAB*mAB) );
  auto eA  = SQRT( pA * pA + m_s1 );
  auto eB  = SQRT( pA * pA + m_s2 );
  auto eC  = ( m_s0 - m_s3 - mAB * mAB ) / ( 2 * mAB );
  auto pC  = SQRT( eC * eC - m_s3 );
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
  /*
  sqCo pos = { gRandom->Uniform(), gRandom->Uniform()};
  double event[12];
  for ( unsigned int i = 0; i < 12; ++i ) event[i] = 0;
  setEvent( pos, event );
  INFO( "s0 = " << m_s0 << " " << sqrt(m_s0) <<  ", x = " << pos.first << " " << pos.second << " r = " << m_min << "  " << m_max );
  double mAB = getMAB(pos);
  double pA  = safe_sqrt( mAB * mAB / 4. - ( m_s1 + m_s2 ) / 2. + ( m_s1 - m_s2 ) * ( m_s1 - m_s2 ) / ( 4 * mAB * mAB ) );
  std::cout << "pA = " << pA << " mAB = " << mAB << std::endl; 
  for ( unsigned int i = 0; i < 12; ++i ) INFO( "Evt[" << i << "] = " << event[i] );
  */
}

/*
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
*/

TGraph* DalitzIntegrator::makeBoundaryGraph(const Projection2D& projection ) const 
{
  AmpGen::Event event(12); 
  auto setEventScalar = [this, &event](const double& f , const double& g  ) mutable 
  {
    sqCo x( f, g ); 
    auto mAB = this->getMAB( x );
    auto pA  = this->safe_sqrt(0.25*mAB*mAB - 0.5*(this->m_s1 + this->m_s2) + 
        (this->m_s1-this->m_s2)*(this->m_s1-this->m_s2)/(4.*mAB*mAB));
    auto eA  = SQRT(pA*pA + m_s1);
    auto eB  = SQRT(pA*pA + m_s2);
    auto eC  = ( this->m_s0 - this->m_s3 - mAB * mAB ) / ( 2. * mAB );
    auto pC  = safe_sqrt( eC * eC - this->m_s3 );
    event[2]  = utils::get<0>(pA);
    event[3]  = utils::get<0>(eA);
    event[6]  = utils::get<0>(-pA);
    event[7]  = utils::get<0>(eB);
    event[9]  = utils::get<0>(pC * sin(M_PI * x.second));
    event[10] = utils::get<0>(pC * cos(M_PI * x.second));
    event[11] = utils::get<0>(eC);
  };
  TGraph* gr = new TGraph(); 
  double eps = 1e-4;
  for( double x = eps; x < 1-eps; x+=0.005)
  {
    setEventScalar( x, eps); 
    auto [px,py] = projection(event); 
    gr->SetPoint( gr->GetN(), px, py); 
  }
  for( double x = eps; x < 1-eps; x+=0.005)
  {
    setEventScalar( 1 -eps,x); 
    auto [px,py] = projection(event); 
    gr->SetPoint( gr->GetN(), px, py); 
  }
  for( double x = eps; x < 1-eps; x+=0.005)
  {
    setEventScalar( 1-x, 1-eps); 
    auto [px,py] = projection(event); 
    gr->SetPoint( gr->GetN(), px, py); 
  }
  for( double x = eps; x < 1-eps; x+=0.005)
  {
    setEventScalar( eps, 1-x); 
    auto [px,py] = projection(event); 
    gr->SetPoint( gr->GetN(), px, py); 
  }
  return gr; 
}

