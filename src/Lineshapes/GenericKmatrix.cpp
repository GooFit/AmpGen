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
namespace AmpGen{ make_enum(PA_TYPE, PVec, PVecPDG, PVecACCMOR, QVec); }

DEFINE_LINESHAPE(GenericKmatrix)
{
  auto props         = ParticlePropertiesList::get( particleName );
  Expression mass    = Parameter( particleName + "_mass", props->mass() );
  Expression radius           = Parameter( particleName + "_radius", props->radius() );
  unsigned nPoles    = NamedParameter<unsigned>(particleName+"::kMatrix::nPoles");
  auto channels      = NamedParameter<std::string>(particleName+"::kMatrix::channels").getVector();
  auto type      = NamedParameter<std::string>(particleName+"::kMatrix::type","default");
  auto kMatrix_formFactor = NamedParameter<std::string>(particleName+"::kMatrix::kMatrix_formFactor","NFF");
  auto pVector_formFactor = NamedParameter<std::string>(particleName+"::kMatrix::pVector_formFactor","NFF");
  auto const pa_type = NamedParameter<PA_TYPE>(particleName+"::kMatrix::production_amplitude",PA_TYPE::PVec);
  auto nChannels     = channels.size();
  auto s0            = mass*mass;
            
  DEBUG( "GenericKmatrix modifier " << lineshapeModifier << " particle = " << particleName );
  auto tokens = split( lineshapeModifier, '.' );
  DEBUG( "GenericKmatrix modifier = " << lineshapeModifier << " nTokens = " << tokens.size() );
  unsigned int cTerm                = tokens.size() > 1 ? stoi( tokens[1] )-1 : 0;
    
  std::vector<Expression> phsps, bw_phase_space, Bls;
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
    Expression gamma  = Parameter(stub + "::gamma");
    DEBUG( "Will link to parameter: " << stub + "::mass");
    poleConfig thisPole(mass*mass,{},gamma,type);
    if( dbexpressions != nullptr ) dbexpressions->emplace_back(stub+"::mass", mass);
    Expression bw_width  = 0;
    Expression bw_width0 = 0;
    for (unsigned channel = 1; channel <= nChannels; ++channel ){
      Expression g = Parameter(stub+"::g::"+std::to_string(channel));
      DEBUG("Will link to parameter: " << stub+"::g::"+std::to_string(channel) );
        
      Expression Bl = 1;
      if(kMatrix_formFactor == "BL"){
            Particle p( channels[channel-1] );
            //Taken from arXiv:1111.6307v1 Eqs. (41), (45)
            Expression k2 = norm(phsps[channel-1]) *s/4.;
            Expression k20 = norm( phaseSpace(mass*mass, p, p.L() )   ) *mass*mass/4.;   // ???
            //Expression k20 = norm( phaseSpace(s0, p, p.L() )   ) * s0/4.;   // ???
            Bl = sqrt( BlattWeisskopf( k2* radius * radius, p.L() ) );
            if(pole==1)Bls.push_back(Bl);
            Bl = Bl / sqrt( BlattWeisskopf( k20 * radius * radius, p.L() ) ) ;
      }
      if(kMatrix_formFactor == "ACCMOR"){
            Particle p( channels[channel-1] );
            //Taken from https://arxiv.org/pdf/0909.2171v1.pdf Eqs. (41), (45)
            Expression k2 = norm(phsps[channel-1]) *s/4.;
            //Expression k20 = norm( phaseSpace(mass*mass, p, p.L() )   ) *s0/4.;
            Bl = sqrt( fpow(k2/(1+k2),p.L()) );
            if(pole==1)Bls.push_back(Bl);
      }
        
      thisPole.add(g, Bl);
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
  
  //Form factor
  const Expression q2         = make_cse( Q2( s, s1, s2 ) );
  const Expression q20        = make_cse( Q2( s0, s1, s2 ) );
  Expression FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( pVector_formFactor == "BL" ) FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  if ( pVector_formFactor == "NFF" ) FormFactor = 1;
  if ( pVector_formFactor == "BELLE2018" ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, L ) );

  Expression I = Constant(0,1);
    
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
        P[k] += (beta * poleConfigs[alpha].couplings[k])/poleConfigs[alpha].pole(s) + phi[k] * (1.+s_0)/(s-s_0);
      }
      F_0 += propagator[{cTerm,k}] * P[k];
    }
    return F_0*FormFactor;
  }
  //P-vector from http://pdg.lbl.gov/2019/reviews/rpp2019-rev-resonances.pdf
  else if(pa_type==PA_TYPE::PVecPDG){
      INFO("Using PDG P-vector approach to build the production amplitude");
      Expression F_0 = 0;

      std::vector<Expression> P(nChannels,0), b(nChannels,0);

      Tensor a( Tensor::dim(nPoles, nChannels) );
      for (unsigned pole = 1; pole <= nPoles; ++pole ){
        for( unsigned channel = 1 ; channel <= nChannels; ++channel ){
            Expression a_Re = Parameter(particleName+"::pole::" + std::to_string(pole)+"::a_re::"+std::to_string(channel));
            Expression a_Im = Parameter(particleName+"::pole::" + std::to_string(pole)+"::a_im::"+std::to_string(channel));
            a[{pole-1,channel-1}] = a_Re + I * a_Im;
        }
      }

      for(unsigned k = 0 ; k < nChannels; ++k){
          Expression b_Re = Parameter(particleName+"::b_re::"+std::to_string(k+1));
          Expression b_Im = Parameter(particleName+"::b_im::"+std::to_string(k+1));
          b[k] = (b_Re + I * b_Im); //  * Bls[k]; ???
      }
      
      for(unsigned k = 0 ; k < nChannels; ++k){
        for(unsigned alpha = 0; alpha < nPoles; ++alpha){
            P[k] += (a[{alpha,k}] * poleConfigs[alpha].couplings[k])/poleConfigs[alpha].pole(s);
        }
        P[k] += b[k];
        F_0 += propagator[{cTerm,k}] * P[k];
      }
      return F_0*FormFactor;
  }
  //P-vector from https://arxiv.org/pdf/0909.2171v1.pdf
  else if(pa_type==PA_TYPE::PVecACCMOR){
      INFO("Using ACCMOR P-vector approach to build the production amplitude");
      Expression F_0 = 0;

      std::vector<Expression> P(nChannels,0), b(nChannels,0);
      Expression tau = Parameter(particleName+"::tau");

      Tensor a( Tensor::dim(nPoles, nChannels) );
      for (unsigned pole = 1; pole <= nPoles; ++pole ){
        for( unsigned channel = 1 ; channel <= nChannels; ++channel ){
            Expression a_Re = Parameter(particleName+"::pole::" + std::to_string(pole)+"::a_re::"+std::to_string(channel));
            Expression a_Im = Parameter(particleName+"::pole::" + std::to_string(pole)+"::a_im::"+std::to_string(channel));
            a[{pole-1,channel-1}] = a_Re + I * a_Im;
        }
      }

      double mKPlus  = 0.493677;
      for(unsigned k = 0 ; k < nChannels; ++k){
          Expression b_Re = Parameter(particleName+"::b_re::"+std::to_string(k+1));
          Expression b_Im = Parameter(particleName+"::b_im::"+std::to_string(k+1));
          b[k] = (b_Re + I * b_Im)/(s-mKPlus*mKPlus);
      }
            
      for(unsigned k = 0 ; k < nChannels; ++k){
        // R_i
        for(unsigned alpha = 0; alpha < nPoles; ++alpha){
            P[k] += (a[{alpha,k}] * poleConfigs[alpha].couplings[k])/poleConfigs[alpha].pole(s);
        }
        // (1 + tau K_ij) D_j
        for(unsigned l = 0 ; l < nChannels; ++l){
            P[k] += ( ( k == l ? 1. : 0. ) + tau * kMatrix[{k, l}] ) * b[l];
        }
        
        F_0 += propagator[{cTerm,k}] * P[k];
      }
      return F_0*FormFactor;
  }
    
  else if(pa_type==PA_TYPE::QVec){
    INFO("Using Q-vector approach to build the production amplitude");
    Expression M;
    for(unsigned k = 0 ; k < nChannels; ++k) M = M + kMatrix[{k,cTerm}] * propagator[{cTerm,k}];
    return M ; // * phsps[0];
  }
  else throw std::runtime_error("This shouldn't happen. Currently supported types: P-vector approach (PA_TYPE::PVec) and Q-vector approach (PA_TYPE::QVec)");
}
