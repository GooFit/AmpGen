#include <stddef.h>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include <exception>

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
#include "AmpGen/enum.h"

using namespace AmpGen;
using namespace AmpGen::fcn;
namespace AmpGen{ make_enum(PA_TYPE, PVec, QVec); }

DEFINE_LINESHAPE(GenericKmatrix)
{
  auto props         = ParticlePropertiesList::get( particleName );
  Expression mass    = Parameter( particleName + "_mass", props->mass() );
  unsigned nPoles    = NamedParameter<unsigned>(particleName+"::kMatrix::nPoles");
  auto channels      = NamedParameter<std::string>(particleName+"::kMatrix::channels").getVector();
  auto const pa_type = NamedParameter<PA_TYPE>(particleName+"::kMatrix::production_amplitude",PA_TYPE::PVec);
  auto nChannels     = channels.size();
  auto s0            = mass*mass;
  std::vector<Expression> phsps, bw_phase_space;
  ADD_DEBUG(s, dbexpressions );
  ADD_DEBUG(s0, dbexpressions );
  INFO("Initialising K-matrix with [nChannels = " << nChannels << ", nPoles = " << nPoles << "]");
  //phase-space
  for( unsigned i = 0 ; i < channels.size(); i+=1 ){
    Particle p( channels[i] );
    INFO( p.decayDescriptor() );
    phsps.emplace_back( phaseSpace(s, p, p.L() ) );
    bw_phase_space.emplace_back( phaseSpace(s0, p, p.L() ) );
    if( dbexpressions != nullptr ) dbexpressions->emplace_back("phsp_"+p.decayDescriptor(), *phsps.rbegin() ); //ADD_DEBUG( *phsps.rbegin(), dbexpressions);
//    ADD_DEBUG( phaseSpace(s0,p,p.L()), dbexpressions );
  }
  //pole configuration for kMatrix (see e.g. eq. (48.25) in http://pdg.lbl.gov/2019/reviews/rpp2019-rev-resonances.pdf)
  std::vector<poleConfig> poleConfigs;
  for (unsigned pole = 1; pole <= nPoles; ++pole ){
    std::string stub = particleName+"::pole::" + std::to_string(pole);
    Expression mass  = Parameter(stub + "::mass");
    DEBUG( "Will link to parameter: " << stub + "::mass");
    poleConfig thisPole(mass*mass);
    if( dbexpressions != nullptr ) dbexpressions->emplace_back(stub+"::mass", mass);
    Expression bw_width  = 0;
    Expression bw_width0 = 0;
    for (unsigned channel = 1; channel <= nChannels; ++channel ){
      Expression g = Parameter(stub+"::g::"+std::to_string(channel));
      DEBUG("Will link to parameter: " << stub+"::g::"+std::to_string(channel) );
      thisPole.add(g, 1);
      if(dbexpressions != nullptr) dbexpressions->emplace_back( stub+"::g::"+std::to_string(channel), g);
      bw_width  = bw_width  + g*g*phsps[channel-1] / mass;
      bw_width0 = bw_width0 + g*g*bw_phase_space[channel-1] / mass;
    }
    if( dbexpressions != nullptr ){
      for( unsigned channel = 1 ; channel <= nChannels; ++channel ){
        Expression g = Parameter(stub+"::g::"+std::to_string(channel));
        Expression BR = g*g*bw_phase_space[channel-1] / ( mass * bw_width0 );
        ADD_DEBUG( BR, dbexpressions );
      }
      ADD_DEBUG(bw_width, dbexpressions);
      ADD_DEBUG(bw_width0, dbexpressions);
    }
    poleConfigs.push_back(thisPole);
  }
  //add non-resonant term to kMatrix (eq. (48.25) in http://pdg.lbl.gov/2019/reviews/rpp2019-rev-resonances.pdf)
  Tensor non_resonant( Tensor::dim(nChannels, nChannels) );
  for(unsigned ch1 = 1; ch1 <= nChannels; ++ch1){
    for( unsigned ch2 = 1; ch2 <= nChannels; ++ch2 ){
      auto c1 = std::to_string(ch1);
      auto c2 = std::to_string(ch2);
      if( ch1 > ch2 ) std::swap(c1,c2);
      std::string nrShape = NamedParameter<std::string>(particleName+"::"+c1+"::"+c2+"::nrShape", "flat");
      Expression f1 = Parameter(particleName+"::f1::"+c1+"::"+c2, 0);
      if( nrShape == "flat") non_resonant[{ch1-1,ch2-1}] = f1;
      else if( nrShape == "pole"){
        Expression f2 = Parameter(particleName+"::f2::"+c1+"::"+c2, 0);
        Expression s0 = Parameter(particleName+"::s0::"+c1+"::"+c2, 0);
        non_resonant[{ch1-1,ch2-1}] = (f1 + f2*sqrt(s)) / (s-s0);
      }
      else WARNING("Unknown shape: " << nrShape);
    }
  }

  Tensor kMatrix = constructKMatrix(s, nChannels, poleConfigs);

  ADD_DEBUG_TENSOR(kMatrix, dbexpressions);
  kMatrix = kMatrix + non_resonant;

  Tensor propagator = getPropagator(kMatrix, phsps);
  ADD_DEBUG_TENSOR(non_resonant, dbexpressions);

  //we have all ingredients to build the production amplitude now
  //follow https://doi.org/10.1007/s1010502a0002 eqns (9)-(13) modulo the Adler zero term (which is argued away in a fuzzy manner by citing a private communication)
  if(pa_type==PA_TYPE::PVec){
    std::vector<Expression> P(nChannels,0), a(nChannels,0), phi(nChannels,0);//the P-vector, a and phi coefficients
    Expression s_0 = Parameter(particleName+"::s0");
    Expression F_0 = 0;//the object we'll return later: the production amplitude in the 0th channel
    //get the coefficients first
    for(unsigned k = 0 ; k < nChannels; ++k) a[k] = Parameter(particleName+"::a::"+std::to_string(k+1));
    //now start loop to calculate the production amplitude
    for(unsigned k = 0 ; k < nChannels; ++k){
      for(unsigned alpha = 0; alpha < nPoles; ++alpha){
        Expression beta = 0;
        for(unsigned q = 0 ; q < nChannels; ++q){
          beta += a[q] * poleConfigs[alpha].couplings[q];
          phi[k] += a[q] * non_resonant[{k,q}];
        }
        P[k] += (beta * poleConfigs[alpha].couplings[k])/(poleConfigs[alpha].s - s) + phi[k] * (1.+s_0)/(s-s_0);
      }
      F_0 += propagator[{0,k}] * P[k];
    }
    //TODO: implement higher orbital angular momentum
    return F_0;
  }
  else if(pa_type==PA_TYPE::QVec){
    Expression M;
    for(unsigned i = 0 ; i < nChannels; ++i) M = M + kMatrix[{i,0}] * propagator[{0,i}];
    return M ; // * phsps[0];
  }
  else throw std::runtime_error("This shouldn't happen. Currently supported types: P-vector approach (PA_TYPE::PVec) and Q-vector approach (PA_TYPE::QVec)");
}
