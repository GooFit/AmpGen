#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Factory.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

struct CoupledChannelConfig {
  std::string p1;
  std::string p2;
  Parameter coupling;  
  size_t L; 
};


DEFINE_LINESHAPE( CoupledChannel )
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  Expression radius = Parameter( particleName + "_radius", props->radius() );

  std::vector< std::string > coupledChannelTokens = NamedParameter<std::string>( particleName + "::coupledChannels" ).getVector();
  INFO( "Got: " << coupledChannelTokens.size() << " coupled channel tokens" );
  if( coupledChannelTokens.size() % 4 != 0 )
  {
    ERROR("Could not parse coupled channel tokens (size = " << coupledChannelTokens.size() << ")" );
    return 0;
  }
  std::vector<CoupledChannelConfig> configs;
  for( size_t x = 0 ; x < coupledChannelTokens.size()/4; ++x ){
    CoupledChannelConfig config;
    config.p1 =  coupledChannelTokens[4*x+0] ;
    config.p2 =  coupledChannelTokens[4*x+1] ;
    config.coupling = Parameter( coupledChannelTokens[4*x+2] );
    config.L  = stoi( coupledChannelTokens[4*x+3] ); 
    configs.push_back( config );
  }
  for( auto& c : configs ){
    INFO( c.p1 << " " << c.p2 << " L = " << c.L << " g = " << c.coupling );
  } 
  return 0;
}


