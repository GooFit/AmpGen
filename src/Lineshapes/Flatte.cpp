
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticlePropertiesList.h"
#include <complex>

using namespace AmpGen;
using namespace AmpGen::fcn;
using namespace std::complex_literals;

Expression aSqrtTerm( const Expression& s, const Expression& m0 )
{
  Expression a2 = 1.0 - ( 4 * m0 * m0 ) / s;
  return Ternary( a2 > 0., sqrt( a2 ), Constant( 0 ) );
}

Expression fSqrtTerm( const Expression& s, const Expression& m0 )
{
  return complex_sqrt( 1.0 - 4*m0*m0/ s );
}
DEFINE_LINESHAPE( Flatte )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression radius = Parameter( particleName + "_radius", props->radius() );
  Expression width  = 0; 
  double mKPlus  = ParticlePropertiesList::get("K+")->mass();
  double mPiPlus = ParticlePropertiesList::get("pi+")->mass();
  double mEta    = ParticlePropertiesList::get("eta0")->mass();
 
  if( particleName == "f(0)(980)0" ){
    double mPi0    = ParticlePropertiesList::get("pi0")->mass();
    double mK0     = ParticlePropertiesList::get("K0")->mass();
    Expression gPi       = Parameter( "Flatte::gPi", 0.165 );
    Expression gK_by_gPi = Parameter( "Flatte::gK_by_gPi", 4.21 );
    Expression Gpipi = (1./3.)*fSqrtTerm(s, mPi0) + (2./3.)*fSqrtTerm(s, mPiPlus);
    Expression GKK   = (1./2.)*fSqrtTerm(s, mK0)  + (1./2.)*fSqrtTerm(s, mKPlus);
    if ( lineshapeModifier == "CutKK" ) 
      GKK = (1./2.) * aSqrtTerm(s, mK0) + (1./2.)*aSqrtTerm( s, mKPlus );
    ADD_DEBUG( Gpipi, dbexpressions );
    ADD_DEBUG( GKK, dbexpressions );
    width = gPi * ( Gpipi + gK_by_gPi * GKK );
  }
  else if( particleName == "a(0)(980)0" ){
    Expression g2pieta         = Parameter("Flatte::g2pieta", 0.175);
    Expression g2KK_by_g2pieta = Parameter("Flatte::g2KK_by_g2pieta",1.20);
    double M_pieta = (mPiPlus + mEta) * (mPiPlus + mEta);
    double P_pieta = (mPiPlus - mEta) * (mPiPlus - mEta);
    Expression Gpieta  = sqrt((1 - P_pieta/s) * (1 - M_pieta/s));
    Expression GKK     = fSqrtTerm(s,mKPlus); 
    ADD_DEBUG( Gpieta, dbexpressions );
    ADD_DEBUG( GKK   , dbexpressions );
    width = g2pieta * ( Gpieta + g2KK_by_g2pieta * GKK );

  }
  auto D = mass*mass - 1i*mass*width; 
  const Expression BW    = 1. / ( D - s ); 
  ADD_DEBUG( mass, dbexpressions );
  ADD_DEBUG( width, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  if( lineshapeModifier.find("kFactor") != std::string::npos ){
    WARNING("k-factor should not be used with the Flatte");
    const Expression kf = kFactor(mass,width,dbexpressions);
    return kf * BW; 
  }
  return BW;
}
