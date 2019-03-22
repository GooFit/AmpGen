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

using namespace AmpGen;

DEFINE_LINESHAPE( PALANO )
{

  const Expression I      = Constant( 0, 1 );
  const Expression sInGeV = s / ( GeV*GeV );
  const double mK         = ParticlePropertiesList::get( "K+" )->mass()  ; 
  const double mPi        = ParticlePropertiesList::get( "pi+" )->mass() ; 
  const double mEtap      = ParticlePropertiesList::get( "eta0" )->mass(); 

  const Expression sNorm     = ( mK * mK + mPi * mPi );
  const Expression s12       = 0.87753; //
  const Expression s32       = 1.0209;
  const Expression I12_adler = sInGeV / sNorm - s12;
  const Expression I32_adler = sInGeV / sNorm - s32;

  std::vector<poleConfig> poleConfigs = {poleConfig( 1.7991, {0.3139, -0.00775} ),
                                         poleConfig( 8.3627, {1.1804, -0.22335} )};
  Expression rho1 = phsp_twoBody( sInGeV, mK, mPi );
  Expression rho2 = phsp_twoBody( sInGeV, mK, mEtap );

  std::vector<Expression> K11_poly = {-0.1553, 0.0909, 0.8618, 0.0629};
  std::vector<Expression> K22_poly = {-0.0036, 0.2590, 1.6950, 2.2300};
  std::vector<Expression> K12_poly = {0.0738, 0.3866, 1.2195, 0.8390};
  double sTop                      = 5.832;
  double sBot                      = 0.360;
  const Expression X = ( 2 * sInGeV - sTop - sBot ) / ( sTop - sBot );
  Tensor kMatrix                   = constructKMatrix( sInGeV, 2, poleConfigs); 
  
  Tensor scattPart( Tensor::dim(2,2));
  scattPart(0,0) = pol(X, K11_poly);
  scattPart(0,1) = pol(X, K12_poly);
  scattPart(1,0) = pol(X, K12_poly);
  scattPart(1,1) = pol(X, K22_poly);
 
  kMatrix = kMatrix + scattPart; 
  const Expression K11  = I12_adler * kMatrix[{0, 0}];
  const Expression K12  = I12_adler * kMatrix[{0, 1}];
  const Expression K22  = I12_adler * kMatrix[{1, 1}];

  const Expression K32 = I32_adler * pol(X, {-0.04046, 0.08143, -0.08849});

  const Expression detK = K11 * K22 - K12 * K12;

  const Expression del = 1 - rho1 * rho2 * detK - I * ( rho1 * K11 + rho2 * K22 );

  const Expression T11 = ( 1. - I * rho2 * K22 );
  const Expression T22 = ( 1. - I * rho1 * K11 );
  const Expression T12 = ( I * rho2 * K12 );

  const Expression T32 = 1 / ( 1 - I * K32 * rho1 );

  ADD_DEBUG( I12_adler, dbexpressions );
  ADD_DEBUG( I32_adler, dbexpressions );
  ADD_DEBUG( K11, dbexpressions );
  ADD_DEBUG( K12, dbexpressions );
  ADD_DEBUG( K22, dbexpressions );
  ADD_DEBUG( K32, dbexpressions );
  ADD_DEBUG( T32, dbexpressions );
  ADD_DEBUG( detK, dbexpressions );
  ADD_DEBUG( del, dbexpressions );
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
