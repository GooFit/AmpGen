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
#include "AmpGen/Particle.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/kMatrix.h"
#include "AmpGen/CoupledChannel.h"

using namespace AmpGen;
using namespace AmpGen::fcn; 

DEFINE_LINESHAPE(GenericKmatrix)
{
  auto props        = ParticlePropertiesList::get( particleName );
  Expression mass   = Parameter( particleName + "_mass", props->mass() );
  INFO( "kMatrix modifier " << lineshapeModifier << " particle = " << particleName );
  auto tokens = split(lineshapeModifier, '.' );
  DEBUG("kMatrix modifier = " << lineshapeModifier << " nTokens = " << tokens.size() );
  size_t nPoles    = NamedParameter<size_t>(     lineshapeModifier + "::kMatrix::nPoles");
  auto channels    = NamedParameter<std::string>(lineshapeModifier + "::kMatrix::channels").getVector();
  size_t nChannels = channels.size();
  std::vector<Expression> phsps;
  std::vector<Expression> bw_phase_space;
  auto s0 = mass*mass;
  ADD_DEBUG(s, dbexpressions );
  ADD_DEBUG(s0, dbexpressions );
  INFO("Initialising K-matrix with : " << nChannels); 
  for( size_t i = 0 ; i < channels.size(); i+=1 ){
    Particle p( channels[i] ); 
    INFO( p.decayDescriptor() );
    Expression sf = Parameter( lineshapeModifier + "::phsp::sf::"+std::to_string(i+1), 1);
    phsps.emplace_back( sf * phaseSpace(s, p, p.orbital() ) );
    bw_phase_space.emplace_back( sf * phaseSpace(s0, p, p.orbital() ) );
    ADD_DEBUG( *phsps.rbegin(), dbexpressions);
    ADD_DEBUG( phaseSpace(s0,p,p.orbital()), dbexpressions );  
  }
  Tensor non_resonant( Tensor::dim(nChannels, nChannels) );
  std::vector<poleConfig> poleConfigs;
  for (size_t pole = 1; pole <= nPoles; ++pole ){
    std::string stub = lineshapeModifier + "::pole::" + std::to_string(pole);
    Expression mass  = Parameter(stub + "::mass");
    poleConfig thisPole(mass*mass);
    if( dbexpressions != nullptr ) dbexpressions->emplace_back(stub+"::mass", mass);
    Expression bw_width  = 0;
    Expression bw_width0 = 0;
    for (size_t channel = 1; channel <= nChannels; ++channel ) 
    {
      Expression g = Parameter(stub+"::g::"+std::to_string(channel));
      thisPole.add(g, 1);
      if( dbexpressions != nullptr ){
        dbexpressions->emplace_back( stub+"::g::"+std::to_string(channel), g); 
      }
      bw_width  = bw_width  + g*g*phsps[channel-1] / mass;
      bw_width0 = bw_width0 + g*g*bw_phase_space[channel-1] / mass;
    }
    for( size_t channel = 1 ; channel <= nChannels; ++channel ){
      Expression g = Parameter(stub+"::g::"+std::to_string(channel));
      Expression BR = g*g*bw_phase_space[channel-1] / ( mass * bw_width0 );
      ADD_DEBUG( BR, dbexpressions );
    }
     ADD_DEBUG(bw_width, dbexpressions);
     ADD_DEBUG(bw_width0, dbexpressions);
    poleConfigs.push_back(thisPole);
  }
  for(size_t ch1 = 1; ch1 <= nChannels; ++ch1){
    for( size_t ch2 = 1; ch2 <= nChannels; ++ch2 ){
      auto c1 = std::to_string(ch1);
      auto c2 = std::to_string(ch2); 
      if( ch1 > ch2 ) std::swap(c1,c2);
      std::string nrShape = NamedParameter<std::string>(lineshapeModifier +"::"+c1+"::"+c2+"::nrShape", "flat");
      Expression f1 = Parameter(lineshapeModifier+"::f1::"+c1+"::"+c2, 0);
      if( nrShape == "flat") non_resonant[{ch1-1,ch2-1}] = f1;
      else if( nrShape == "pole"){ 
        Expression f2 = Parameter(lineshapeModifier+"::f2::"+c1+"::"+c2, 0);
        Expression s0 = Parameter(lineshapeModifier+"::s0::"+c1+"::"+c2, 0);
        non_resonant[{ch1-1,ch2-1}] = (f1 + f2*sqrt(s)) / (s-s0);
      }
      else WARNING("Unknown shape: " << nrShape); 
    }
  }
  Tensor kMatrix = constructKMatrix(s, nChannels, poleConfigs);

  ADD_DEBUG_TENSOR(kMatrix   , dbexpressions);
  kMatrix = kMatrix + non_resonant; 

  Tensor propagator = getPropagator(kMatrix, phsps);
  ADD_DEBUG_TENSOR(non_resonant, dbexpressions);
  Expression M;
  for(size_t i = 0 ; i < nChannels; ++i) M = M + kMatrix[{i,0}] * propagator[{0,i}];
  return M ; // * phsps[0];
}
