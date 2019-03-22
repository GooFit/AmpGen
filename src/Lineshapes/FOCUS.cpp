#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/kMatrix.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

DEFINE_LINESHAPE( FOCUS )
{
  const Expression I      = Constant( 0, 1 );
  const Expression sInGeV = s / (GeV*GeV);
  const double mK         = ParticlePropertiesList::get( "K+" )->mass() ;
  const double mPi        = ParticlePropertiesList::get( "pi+" )->mass() ;
  const double mEtap      = ParticlePropertiesList::get( "eta'(958)0" )->mass();

  const Expression sNorm     = ( mK * mK + mPi * mPi );
  const Expression s12       = 0.23;
  const Expression s32       = 0.27;
  const Expression I12_adler = ( sInGeV - s12 ) / sNorm;
  const Expression I32_adler = ( sInGeV - s32 ) / sNorm;

  std::vector<poleConfig> poleConfigs = {poleConfig( 1.7919, {0.31072, -0.02323} )};
  Expression rho1                     = phsp_FOCUS( sInGeV, mK, mPi );
  Expression rho2                     = phsp_FOCUS( sInGeV, mK, mEtap );
  const Expression X = ( sInGeV / sNorm ) - 1;

  Tensor kMatrix = constructKMatrix( sInGeV, 2, poleConfigs); 
  Tensor scattPart( Tensor::dim(2,2));
  scattPart(0,0) = pol(X, {0.79299, -0.15099, 0.00811} );
  scattPart(1,1) = pol(X, {0.17054, -0.0219, 0.00085655} );
  scattPart(0,1) = pol(X, {0.15040, -0.038266, 0.0022596} );
  scattPart(1,0) = pol(X, {0.15040, -0.038266, 0.0022596} );

  kMatrix = kMatrix + scattPart;   
  const Expression K11 = I12_adler * kMatrix[{0, 0}];
  const Expression K12 = I12_adler * kMatrix[{0, 1}];
  const Expression K22 = I12_adler * kMatrix[{1, 1}];

  const Expression K32 = I32_adler * pol( X, {-0.22147, 0.026637, -0.00092057} );

  const Expression detK = K11 * K22 - K12 * K12;

  const Expression del = 1 - rho1 * rho2 * detK - I * ( rho1 * K11 + rho2 * K22 );

  const Expression T11 = ( 1. - I * rho2 * K22 );
  const Expression T22 = ( 1. - I * rho1 * K11 );
  const Expression T12 = ( I * rho2 * K12 );

  const Expression T32 = 1 / ( 1 - I * K32 * rho1 );
  ADD_DEBUG(sInGeV   , dbexpressions );
  ADD_DEBUG(mK       , dbexpressions ); 
  ADD_DEBUG(sNorm    , dbexpressions );
  ADD_DEBUG(X        , dbexpressions );
  ADD_DEBUG(I12_adler, dbexpressions );
  ADD_DEBUG(I32_adler, dbexpressions );
  ADD_DEBUG(K11      , dbexpressions );
  ADD_DEBUG(K12      , dbexpressions );
  ADD_DEBUG(K22      , dbexpressions );
  ADD_DEBUG(K32      , dbexpressions );
  ADD_DEBUG(T32      , dbexpressions );
  ADD_DEBUG(detK     , dbexpressions );
  ADD_DEBUG(del      , dbexpressions );
  if ( lineshapeModifier == "Kpi" )
    return ( K11 - I * rho2 * detK ) / del;
  else if ( lineshapeModifier == "KEta" )
    return K12 / del;
  else if ( lineshapeModifier == "I32" )
    return T32;
  else {
    ERROR( "P-vector component : " << lineshapeModifier << " is not recognised" );
    return 1;
  }
}
