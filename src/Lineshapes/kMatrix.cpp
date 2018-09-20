#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/kMatrix.h"
#include "AmpGen/Units.h"

using namespace AmpGen;

Expression AmpGen::phsp_twoBody( const Expression& s, const double& m0, const double& m1 )
{
  return fcn::complex_sqrt( 1.0 - ( m0 + m1 ) * ( m0 + m1 ) / s );
}

Expression AmpGen::phsp_fourPi( const Expression& s )
{
  double mPiPlus( 0.139570 );
  Expression rho_4pi = pol( s, {0.00051, -0.01933, 0.13851, -0.20840, -0.29744, 0.13655, 1.07885} );
  return Ternary( s > 1, phsp_twoBody( s, 2 * mPiPlus, 2 * mPiPlus ), rho_4pi );
}

Expression AmpGen::phsp_FOCUS( const Expression& s, const double& m0, const double& m1 )
{
  double mp  = ( m0 + m1 );
  double mm  = ( m0 - m1 );
  return fcn::complex_sqrt( ( 1.0 - mp * mp / s ) * ( 1.0 - mm * mm / s ) ); 
}

std::vector<Parameter> AmpGen::paramVector( const std::string& name, const unsigned int& nParam )
{
  std::vector<Parameter> returnVector;
  for ( unsigned int i = 0; i < nParam; ++i ) returnVector.emplace_back( name + std::to_string( i ) );
  return returnVector;
}

Expression AmpGen::gFromGamma( const Expression& m, const Expression& gamma, const Expression& rho )
{
  return Sqrt( m * gamma / rho );
}

Tensor AmpGen::getPropagator( const Tensor& kMatrix, const std::vector<Expression>& phaseSpace,
                              const Expression& adlerTerm )
{
  unsigned int nChannels = kMatrix.dims()[0];
  Tensor TMatrx( std::vector<size_t>( {nChannels, nChannels} ) );
  DEBUG( "Constructing T-Matrix" );
  Expression J = Constant(0,1);
  for ( unsigned int i = 0; i < nChannels; ++i ) {
    for ( unsigned int j = 0; j < nChannels; ++j ) {
      TMatrx[{i, j}] = SubTree( ( i == j ? 1 : 0 ) - J * adlerTerm  * kMatrix[{i, j}] * phaseSpace[j] );
    }
  }
  return TMatrx.Invert();
}

Tensor AmpGen::constructKMatrix(
    const Expression& this_s, const unsigned int& nChannels, const std::vector<poleConfig>& poleConfigs,
    const std::function<Expression( const unsigned int& i, const unsigned int& j, const Expression& s )>& SVP,
    DebugSymbols* dbexpressions )
{

  Tensor kMatrix( std::vector<size_t>( {nChannels, nChannels} ) );
  for ( unsigned int i = 0; i < nChannels; ++i ) {
    for ( unsigned int j = 0; j < nChannels; ++j ) {
      Expression sumOverPoles = 0;
      for ( auto& pole : poleConfigs ) {
        Expression term = ( pole.couplings[i] * pole.couplings[j] ) / ( pole.s - this_s );
        ADD_DEBUG( term, dbexpressions );
        sumOverPoles = sumOverPoles + term;
      }
      Expression nr = SVP == nullptr ? Expression( 0 ) : SVP( i, j, this_s );
      if ( dbexpressions != nullptr )
        dbexpressions->emplace_back( "SVP[" + std::to_string( i ) + "," + std::to_string( j ) + "]", nr );
      kMatrix[{i, j}] = sumOverPoles + nr;
    }
  }
  return kMatrix;
}

DEFINE_LINESHAPE( kMatrix )
{
  Expression sInGeV = SubTree(s) / (GeV*GeV);
  DEBUG( "kMatrix modifier " << lineshapeModifier << " particle = " << particleName );
  auto tokens = split( lineshapeModifier, '.' );

  DEBUG( "kMatrix modifier = " << lineshapeModifier << " nTokens = " << tokens.size() );
  unsigned int pTerm                = tokens.size() > 1 ? stoi( tokens[1] ) : 999;
  unsigned int nPoles               = 5;
  unsigned int nChannels            = 5;
  std::vector<std::string> channels = {"pipi", "KK", "4pi", "EtaEta", "EtapEta"};

  double mPiPlus = 0.139570;
  double mKPlus  = 0.493677;
  double mEta    = 0.547862;
  double mEtap   = 0.96778;

  Expression sA0      = Parameter( "sA0", -0.15 );
  Expression sA       = Parameter( "sA", 1.0 );
  Expression s0_prod  = Parameter( particleName + "_s0_prod", -0.07 );
  Expression s0_scatt = Parameter( "s0_scatt", -3.92637 );

  std::vector<Parameter> fScatt = paramVector( "f_scatt", nChannels );
  std::vector<poleConfig> poleConfigs; 

  for ( unsigned int pole = 1; pole <= nPoles; ++pole ) {
    std::string stub = "IS_p" + std::to_string( pole ) + "_";
    Expression mass  = Parameter( stub + "mass" );
    poleConfig p( mass * mass );
    for ( unsigned int ch = 0; ch < nChannels; ++ch ) p.add( Parameter( stub + channels[ch] ) );
    poleConfigs.push_back( p );
  }

  std::vector<Expression> phaseSpace = {phsp_twoBody( sInGeV, mPiPlus, mPiPlus ),
                                        phsp_twoBody( sInGeV, mKPlus, mKPlus ), phsp_fourPi( sInGeV ),
                                        phsp_twoBody( sInGeV, mEta, mEta ), phsp_twoBody( sInGeV, mEta, mEtap )};

  auto scatteringTerm = [&]( const unsigned int& i, const unsigned int& j, const Expression& s ) {
    if ( i == 0 ) return fScatt[j] * ( 1 - s0_scatt ) / ( s - s0_scatt );
    if ( j == 0 ) return fScatt[i] * ( 1 - s0_scatt ) / ( s - s0_scatt );
    return Expression( 0 );
  };

  auto kMatrix = constructKMatrix( sInGeV, 5, poleConfigs, scatteringTerm );

  kMatrix.imposeSymmetry(0,1);

  Expression adlerTerm = SubTree( ( 1. - sA0 ) * ( sInGeV - sA * mPiPlus * mPiPlus / 2 ) / ( sInGeV - sA0 ) );
  Tensor F             = getPropagator( kMatrix, phaseSpace, adlerTerm );

  if ( dbexpressions != nullptr ) {
    for ( auto& p : poleConfigs ) {
      ADD_DEBUG( p.s, dbexpressions );
      ADD_DEBUG( p.couplings[0], dbexpressions );
      ADD_DEBUG( p.couplings[1], dbexpressions );
      ADD_DEBUG( p.couplings[2], dbexpressions );
      ADD_DEBUG( p.couplings[3], dbexpressions );
      ADD_DEBUG( p.couplings[4], dbexpressions );
    }
    for ( unsigned int i = 0; i < 5; ++i )
      dbexpressions->emplace_back( "phi[" + std::to_string( i ) + "]", phaseSpace[i] );

    ADD_DEBUG( sA * mPiPlus * mPiPlus / 2, dbexpressions );
    ADD_DEBUG( sA0, dbexpressions );
    ADD_DEBUG( adlerTerm, dbexpressions );
    for ( unsigned int i = 0; i < 5; ++i ) {
      for ( unsigned int j = 0; j < 5; ++j ) {
        dbexpressions->emplace_back( "K[" + std::to_string( i ) + "," + std::to_string( j ) + "]", kMatrix[{i, j}] );
      }
    }
  }

  if ( tokens[0] == "scatt" ) {
    unsigned int i = stoi( tokens[1] );
    unsigned int j = stoi( tokens[2] );
    Expression M;
    for ( unsigned int x = 0; x < nChannels; ++x ) {
      M = M + kMatrix[{x, i}] * F[{j, x}];
    }
    return M;
  } else if ( tokens[0] == "pole" ) {
    Expression M;
    DEBUG( "Returning pole term" );
    auto pole = poleConfigs[pTerm];
    for ( unsigned int i = 0; i < pole.couplings.size(); ++i ) {
      M = M + F[{0, i}] * pole.couplings[i];
    }
    return SubTree( M / ( pole.s - sInGeV ) );
  } else if ( tokens[0] == "poleKK" ) {
    Expression M;
    DEBUG( "Returning pole term" );
    auto pole = poleConfigs[pTerm];
    for ( unsigned int i = 0; i < pole.couplings.size(); ++i ) {
      M = M + F[{1, i}] * pole.couplings[i];
    }
    return SubTree( M / ( pole.s - sInGeV ) );
  }

  else if ( tokens[0] == "prod" ) {
    DEBUG( "Returning bkg term" );
    return F[{0, pTerm}] * ( 1 - s0_prod ) / ( sInGeV - s0_prod );
  } else if ( tokens[0] == "prodKK" ) {
    DEBUG( "Returning bkg term" );
    return F[{1, pTerm}] * ( 1 - s0_prod ) / ( sInGeV - s0_prod );
  }
  ERROR( "Lineshape not found: " << lineshapeModifier );
  return Expression( 0 );
}
