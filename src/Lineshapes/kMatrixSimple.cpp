#include <stddef.h>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/kMatrix.h"

using namespace AmpGen;

DEFINE_LINESHAPE( kMatrixSimple )
{
  Expression sInGeV = SubTree( s/(GeV*GeV) );
  DEBUG( "kMatrix modifier " << lineshapeModifier << " particle = " << particleName );
  auto tokens = split( lineshapeModifier, '.' );

  DEBUG( "kMatrix modifier = " << lineshapeModifier << " nTokens = " << tokens.size() );
  unsigned int pTerm = tokens.size() > 0 ? stoi( tokens[0] ) : 999;
  if ( pTerm == 999 ) {
    ERROR( "Pole term not recognised: " << pTerm );
  }

  size_t nPoles    = NamedParameter<size_t>( "kMatrix::nPoles", 2 );
  size_t nChannels = NamedParameter<size_t>( "kMatrix::nChannels", 2 );

  std::vector<Expression> phaseSpace;
  std::vector<std::pair<double, double>> masses;
  for ( unsigned int i = 1; i <= nChannels; ++i ) {
    std::vector<std::string> particlesInThisChannel = NamedParameter<std::string>( "kMatrix::Channel::" + std::to_string(i) ).getVector();
    if ( particlesInThisChannel.size() != 2 ) ERROR( "Only does two body channels for now" );
    double m1 = ParticlePropertiesList::get( particlesInThisChannel[0] )->mass();
    double m2 = ParticlePropertiesList::get( particlesInThisChannel[1] )->mass();
    masses.emplace_back( m1, m2 );
    phaseSpace.push_back( phsp_twoBody( sInGeV, m1, m2 ) );
  }

  std::vector<poleConfig> poleConfigs;

  for ( unsigned int pole = 1; pole <= nPoles; ++pole ) {
    std::string stub = "kMatrix::pole::" + std::to_string( pole );
    Expression mass  = Parameter( ( stub + "::mass" ) );
    poleConfig thisPole( mass * mass );
    for ( unsigned int channel = 0; channel < nChannels; ++channel ) {
      Expression gPiPi = Parameter( ( stub + "::width::" + std::to_string( channel ) ) );
      Expression rho0  = phsp_twoBody( mass * mass, masses[channel].first, masses[channel].second );
      Expression g     = gFromGamma( mass, gPiPi, rho0 );
      thisPole.add( g );
      if ( dbexpressions != nullptr ) {
        Expression rhoRatio = phaseSpace[channel] / rho0;
        Expression GammaS = g * g * phaseSpace[channel] / mass;
        ADD_DEBUG( gPiPi, dbexpressions );
        ADD_DEBUG( g, dbexpressions );
        ADD_DEBUG( rhoRatio, dbexpressions );
        ADD_DEBUG( rho0, dbexpressions );
        ADD_DEBUG( GammaS, dbexpressions );
      }
    }
    poleConfigs.push_back( thisPole );
  }

  auto kMatrix = constructKMatrix( sInGeV, nChannels, poleConfigs);
  ADD_DEBUG_TENSOR( kMatrix, dbexpressions );
  Tensor propagator = getPropagator( kMatrix, phaseSpace );
  Expression M;

  auto pole = poleConfigs[pTerm];
  if ( dbexpressions != nullptr ) dbexpressions->emplace_back( "P[0,0]", propagator[{0, 0}] );
  for ( unsigned int i = 0; i < pole.couplings.size(); ++i ) 
  {
    M = M + propagator[{0, i}] * pole.g(i);
  }
  return SubTree( M / ( pole.s - sInGeV ) );
}
