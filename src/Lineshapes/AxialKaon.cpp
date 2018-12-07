
#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Lineshapes.h"
#include "AmpGen/kMatrix.h"
#include "AmpGen/Spline.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

DEFINE_LINESHAPE( AxialKaon )
{
  ERROR("This lineshape does not work as intended, should not be used");
/*
  
  Spline rho_Krho("rho_Krho",20, 0.597247, 5, 
      {5.35636e-23, 0.000210696, 0.00208479, 0.00993997, 0.0415061, 0.120885, 0.196599, 0.256062, 0.304949, 0.34693, 0.384132, 0.417868, 0.448995, 0.478101, 0.505606, 0.531815, 0.556963, 0.581226, 0.604745, 0.627631});
  
  Spline rho_Kstpi("rho_Kspi",20,0.597247, 5, 
      {6.29163e-23, 0.000390208, 0.0161092, 0.0855491, 0.128145, 0.159565, 0.184867, 0.20627, 0.224985, 0.241753, 0.257056, 0.271229, 0.284507, 0.297066, 0.309035, 0.320516, 0.331588, 0.342313, 0.35274, 0.362912} );    

  Spline rho_KstpiD("rho_KspiD",20,0.597247, 5, 
    {6.29163e-23, 0.000390208, 0.0161092, 0.0855491, 0.128145, 0.159565, 0.184867, 0.20627, 0.224985, 0.241753, 0.257056, 0.271229, 0.284507, 0.297066, 0.309035, 0.320516, 0.331588, 0.342313, 0.35274, 0.362912} );

  Spline rho_KSpi("rho_KSpi",20,0.597247, 5, 
    {5.58343e-23, 0.000178087, 0.00132915, 0.00423057, 0.00963868, 0.0185276, 0.0324146, 0.053891, 0.0870863, 0.13619, 0.202066, 0.282239, 0.373053, 0.470919, 0.572771, 0.67618, 0.779316, 0.880855, 0.979877, 1.07577});


  std::vector<Expression> phaseSpace = { rho_Krho(s), rho_Kstpi(s), rho_KstpiD(s), rho_KSpi(s) };

  std::vector< std::string > channels = {"Krho","Kspi","KspiD","KSpi"};
  unsigned int nPoles = 2;
  std::vector<poleConfig>  poleConfigs;
  for( unsigned int pole=1;pole <= nPoles; ++pole )
  {
    std::string stub = "AxialKaon_p"+std::to_string(pole)+"_";
    Expression mass = Parameter( stub + "mass" );
    poleConfig p( mass * mass );
    for ( unsigned int ch = 0; ch < 4; ++ch ) p.add( Parameter( stub + channels[ch] ) );
    poleConfigs.push_back( p );
  } 
  auto kMatrix = constructKMatrix( s, 4, poleConfigs );

  Tensor F = getPropagator( kMatrix, phaseSpace );

  auto tokens = split( lineshapeModifier, '.' );
  unsigned int pTerm   = stoi(tokens[0]);
  unsigned int channel = stoi(tokens[1]);
  auto pole = poleConfigs[pTerm];
  Expression M = 0;
  for ( unsigned int i = 0; i < pole.couplings.size(); ++i ) {
    M = M + F[{channel, i}] * pole.couplings[i];
  }
  return SubTree( M / ( pole.s - s ) );
  */
  return 0;
}
