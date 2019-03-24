#include "AmpGen/Kinematics.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>
#include <cmath>

#include "AmpGen/Expression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Event.h"

using namespace AmpGen;

MomentumTransfer::MomentumTransfer( const std::vector<size_t>& _p1, const std::vector<size_t>& _p2 )
    : p1( _p1 ), p2( _p2 )
{
  for ( auto& p : p1 ) s.push_back( p );
  for ( auto& p : p2 ) s.push_back( p );
}

double MomentumTransfer::operator()( const Event& evt ) const
{
  double s0 = evt.s( s );
  double s1 = evt.s( p1 );
  double s2 = evt.s( p2 );
  return sqrt( MomentumTransfer::Q2( s0, s1, s2 ) );
}
double MomentumTransfer::Q2( const double& s, const double& s1, const double& s2 ) const
{
  return s / 4. - ( s1 + s2 ) / 2. + ( s1 - s2 ) * ( s1 - s2 ) / ( 4 * s );
}

HelicityCosine::HelicityCosine( const std::vector<size_t>& p1, const std::vector<size_t>& p2,
                                const std::vector<size_t>& pR )
    : _i( p1 ), _j( p2 ), _pR( pR )
{
}

HelicityCosine::HelicityCosine( const size_t& i, const size_t& j, const std::vector<size_t>& pR )
    : _i( 1, i ), _j( 1, j ), _pR( pR )
{
}

double HelicityCosine::operator()( std::vector<Event>::iterator evt ) const { return ( *this )( *evt ); }
double HelicityCosine::operator()( const Event& evt ) const
{
  TLorentzVector PR = pFromEvent( evt, _pR );
  TLorentzVector pi = pFromEvent( evt, _i );
  TLorentzVector pj = pFromEvent( evt, _j );
  return dotProduct(pi, pj, PR) / sqrt( dotProduct( pi, pi, PR ) * dotProduct( pj, pj, PR ) );
}

TLorentzVector AmpGen::pFromEvent( const Event& evt, const size_t& ref )
{
  return TLorentzVector( evt.address( 4 * ref ) );
}

TLorentzVector AmpGen::pFromEvent( const Event& evt, const std::vector<size_t>& ref )
{
  double px( 0 ), py( 0 ), pz( 0 ), pE( 0 );
  for ( auto& r : ref ) {
    px += evt[4 * r + 0];
    py += evt[4 * r + 1];
    pz += evt[4 * r + 2];
    pE += evt[4 * r + 3];
  }
  return TLorentzVector( px, py, pz, pE );
}

double AmpGen::acoplanarity( const Event& evt )
{
  TLorentzVector p0 = pFromEvent( evt, 0 );
  TLorentzVector p1 = pFromEvent( evt, 1 );
  TLorentzVector p2 = pFromEvent( evt, 2 );
  TLorentzVector p3 = pFromEvent( evt, 3 );
  TLorentzVector pD = p0 + p1 + p2 + p3;
  p0.Boost( -pD.BoostVector() );
  p1.Boost( -pD.BoostVector() );
  p2.Boost( -pD.BoostVector() );
  p3.Boost( -pD.BoostVector() );
  TVector3 t1 = ( p0.Vect().Cross( p1.Vect() ) ).Unit();
  TVector3 t2 = ( p2.Vect().Cross( p3.Vect() ) ).Unit();
  return acos( t1.Dot( t2 ) );
}

double AmpGen::dotProduct( const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& pX )
{
  return -p1.Dot( p2 ) + p1.Dot( pX ) * p2.Dot( pX ) / pX.Dot( pX );
}

double AmpGen::PHI( const Event& evt )
{
  TLorentzVector pA_4vec = pFromEvent( evt, 0 );
  TLorentzVector pB_4vec = pFromEvent( evt, 1 );
  TLorentzVector pC_4vec = pFromEvent( evt, 3 );
  TLorentzVector pD_4vec = pFromEvent( evt, 2 );
  TLorentzVector pM      = pA_4vec + pB_4vec + pC_4vec + pD_4vec;
  pA_4vec.Boost( -pM.BoostVector() );
  pB_4vec.Boost( -pM.BoostVector() );
  pC_4vec.Boost( -pM.BoostVector() );
  pD_4vec.Boost( -pM.BoostVector() );
  TVector3 pA = pA_4vec.Vect();
  TVector3 pB = pB_4vec.Vect();
  TVector3 pC = pC_4vec.Vect();
  TVector3 pD = pD_4vec.Vect();
  TVector3 n1 = ( pA.Cross( pB ) ).Unit();
  TVector3 n2 = ( pC.Cross( pD ) ).Unit();

  double cosPhi = n1.Dot( n2 );
  double sinPhi = pD.Dot( n1 ) / ( pD.Cross( ( pA + pB ).Unit() ).Mag() );
  double phi    = TMath::ATan2( sinPhi, cosPhi );
  return phi;
}

double AmpGen::phi( const Event& evt, int i, int j, int k, int w )
{
  TLorentzVector pA_4vec = pFromEvent( evt, i );
  TLorentzVector pB_4vec = pFromEvent( evt, j );
  TLorentzVector pC_4vec = pFromEvent( evt, k );
  TLorentzVector pD_4vec = pFromEvent( evt, w );
  TLorentzVector pM      = pA_4vec + pB_4vec + pC_4vec + pD_4vec;
  pA_4vec.Boost( -pM.BoostVector() );
  pB_4vec.Boost( -pM.BoostVector() );
  pC_4vec.Boost( -pM.BoostVector() );
  pD_4vec.Boost( -pM.BoostVector() );
  TVector3 pA = pA_4vec.Vect();
  TVector3 pB = pB_4vec.Vect();
  TVector3 pC = pC_4vec.Vect();
  TVector3 pD = pD_4vec.Vect();
  TVector3 n1 = ( pA.Cross( pB ) ).Unit();
  TVector3 n2 = ( pC.Cross( pD ) ).Unit();

  double cosPhi = n1.Dot( n2 );
  double sinPhi = pD.Dot( n1 ) / ( pD.Cross( ( pA + pB ).Unit() ).Mag() );
  double phi    = TMath::ATan2( sinPhi, cosPhi );
  return phi;
}

void AmpGen::boost( Event& evt, const std::tuple<double, double, double>& n, const double& v )
{
  double gamma = 1. / sqrt( 1 - v * v );
  double nx   = std::get<0>( n );
  double ny   = std::get<1>( n );
  double nz   = std::get<2>( n );
  double norm = sqrt( nx * nx + ny * ny + nz * nz );

  for ( size_t i = 0; i < evt.size() / 4; ++i ) {
    double nv = evt[4 * i] * nx + evt[4 * i + 1] * ny + evt[4 * i + 2] * nz;
    evt[4*i+0] += ( (gamma-1) * nv / norm + gamma * evt[4*i+3] * v ) * nx / norm;
    evt[4*i+1] += ( (gamma-1) * nv / norm + gamma * evt[4*i+3] * v ) * ny / norm;
    evt[4*i+2] += ( (gamma-1) * nv / norm + gamma * evt[4*i+3] * v ) * nz / norm;
    evt[4*i+3] = gamma * ( evt[4*i+3] + v * nv / norm );
  }
}

void AmpGen::rotate( Event& evt, const std::tuple<double,double,double>& n, const double& v )
{
  double nx   = std::get<0>( n );
  double ny   = std::get<1>( n );
  double nz   = std::get<2>( n );
  double cv   = cos(v);
  double sv   = sin(v);
  double norm = sqrt( nx * nx + ny * ny + nz * nz );
  for( size_t i = 0 ; i < evt.size()/4;++i){
    double ix = evt[ 4*i+0];
    double iy = evt[ 4*i+1];
    double iz = evt[ 4*i+2];
    double k = (1-cv) * (ix*nx + iy*ny +iz*nz) ; 
    evt[4*i + 0] = cv * ix + sv*( iz * ny - iy * nz )/norm + k * nx / ( norm * norm );
    evt[4*i + 1] = cv * iy + sv*( ix * nz - iz * nx )/norm + k * ny / ( norm * norm );
    evt[4*i + 2] = cv * iz + sv*( iy * nx - ix * ny )/norm + k * nz / ( norm * norm );
  }
}

void AmpGen::rotateBasis( Event& evt, const TVector3& p1, const TVector3& p2, const TVector3& p3 )
{
  for(size_t i = 0 ; i < evt.size()/4;++i)
  { 
    double ex = evt[4*i+0];
    double ey = evt[4*i+1];
    double ez = evt[4*i+2];
    evt[4*i+0] = (p1.x() * ex + p1.y() * ey + p1.z() * ez )/p1.Mag();
    evt[4*i+1] = (p2.x() * ex + p2.y() * ey + p2.z() * ez )/p2.Mag();
    evt[4*i+2] = (p3.x() * ex + p3.y() * ey + p3.z() * ez )/p3.Mag();
  }
}
