#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/kMatrix.h"

DEFINE_LINESHAPE( ObelixRho )
{
  Expression sInGeV = SubTree( s / (GeV*GeV) );
  DEBUG( "kMatrix modifier " << lineshapeModifier << " particle = " << particleName );
  auto tokens = split( lineshapeModifier, '.' );

  DEBUG( "kMatrix modifier = " << lineshapeModifier << " nTokens = " << tokens.size() );
  unsigned int pTerm = tokens.size() > 0 ? stoi( tokens[0] ) : 999;
  if ( pTerm == 999 ) {
    ERROR( "Pole term not recognised: " << pTerm );
  }

  double mPiPlus( 0.139570 );
  double mKPlus( 0.493677 );
  std::vector<Expression> phaseSpace = {phsp_twoBody( sInGeV, mPiPlus, mPiPlus ), phsp_fourPi( sInGeV ),
                                        phsp_twoBody( sInGeV, mKPlus, mKPlus )};
  const unsigned int nPoles    = 3;
  const unsigned int nChannels = 3;
  auto props                   = ParticlePropertiesList::get( "rho(770)0" );

  const Expression radius = Parameter( particleName + "_radius", props->radius() );

  std::vector<poleConfig> poleConfigs;

  for ( unsigned int i = 1; i <= nPoles; ++i ) {
    std::string stub = "rho_p" + std::to_string( i );
    Expression mRho  = Parameter( ( stub + "_mass" ) );
    Expression gRho1 = Parameter( ( stub + "_pipi" ) );
    Expression gRho2 = Parameter( ( stub + "_4pi" ) );
    Expression gRho3 = Parameter( ( stub + "_KK" ) );
    poleConfig pole( mRho * mRho );

    pole.add( gFromGamma( mRho, gRho1, phsp_twoBody( mRho * mRho, mPiPlus, mPiPlus ) ),
              BL( sInGeV, mRho * mRho, mPiPlus * mPiPlus, mPiPlus * mPiPlus, radius, 1 ) );

    pole.add( gFromGamma( mRho, gRho2, phsp_fourPi( mRho * mRho ) ) );

    pole.add( gFromGamma( mRho, gRho3, phsp_twoBody( mRho * mRho, mKPlus, mKPlus ) ),
              BL( sInGeV, mRho * mRho, mKPlus * mKPlus, mKPlus * mKPlus, radius, 1 ) );

    poleConfigs.emplace_back( pole );
  }

  auto kMatrix = constructKMatrix( sInGeV, nChannels, poleConfigs );

  if ( dbexpressions != nullptr ) {
    for ( auto& p : poleConfigs ) {
      ADD_DEBUG( p.s, dbexpressions );
      ADD_DEBUG( p.couplings[0], dbexpressions );
      ADD_DEBUG( p.couplings[1], dbexpressions );
      ADD_DEBUG( p.couplings[2], dbexpressions );
    }
    ADD_DEBUG( phaseSpace[0], dbexpressions );
    ADD_DEBUG( phaseSpace[1], dbexpressions );
    ADD_DEBUG( phaseSpace[2], dbexpressions );
    for ( unsigned int i = 0; i < nChannels; ++i ) {
      for ( unsigned int j = 0; j < nChannels; ++j ) {
        dbexpressions->emplace_back( "K[" + std::to_string( i ) + "," + std::to_string( j ) + "]", kMatrix[{i, j}] );
      }
    }
  }

  Tensor propagator = getPropagator( kMatrix, phaseSpace );

  Expression M;

  auto pole           = poleConfigs[pTerm];
  const Expression q2 = Abs( Q2( s, s1, s2 ) );
  //  Expression FormFactor = Sqrt(BlattWeisskopf_Norm( q2*radius*radius,0, L));

  for ( unsigned int i = 0; i < nChannels; ++i ) {
    M = M + propagator[{0, i}] * pole.g( i );
  }
  return SubTree( M / ( pole.s - sInGeV ) ); // /Sqrt(q2) ;
}
