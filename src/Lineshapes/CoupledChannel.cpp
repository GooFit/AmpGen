#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Particle.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/CoupledChannel.h"

using namespace AmpGen;
using namespace AmpGen::fcn; 

Expression H(const Expression& x, 
             const Expression& y, 
             const Expression& z )
{
  auto k = ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*z*y) ; 
  return Ternary( Imag(sqrt(k)) > 0 , sqrt(k) , -sqrt(k) ); /// change the branch of the sqrt
}

Expression Hr(const Expression& x, 
              const Expression& y, 
              const Expression& z )
{
  auto k = ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*z*y) ; 
  return complex_sqrt(k);
}

Expression rho_twoBody( const Expression& s, const Expression& s1, const Expression& s2)
{
  return Hr(s,s1,s2)/(16*M_PI*s);
}

Expression rho_threeBody( const Expression& s, 
                          const Particle& resonance, 
                          const Particle& bachelor ) 
{
  const ParticleProperties& is = *resonance.props();
  const ParticleProperties& p1 = *bachelor.props();
  const ParticleProperties& p2 = *resonance.daughter(0)->props();
  const ParticleProperties& p3 = *resonance.daughter(1)->props();
    
  Expression s1 = p1.mass() * p1.mass();
  Expression s2 = p2.mass() * p2.mass();
  Expression s3 = p3.mass() * p3.mass();
  Expression sT =  s2 + s3 + 2*sqrt(s2*s3); 
  Expression s23_max =  s  + s1 - 2*sqrt(s*s1);
  Expression i = Constant(0,1);
  Expression m = is.mass();
  Expression G = is.width();
  auto z   = is.mass() * ( is.mass() - i * is.width() ); 
  auto Hz  = H(s,z,s1);
  auto HsT = Hr(s,sT,s1);
  auto K   = Hz*HsT + z*(sT-s-s1) + (s-s1)*(s-s1) - sT*(s+s1);

  return    Ternary( s23_max > sT, Imag( 
      z * log( (s+s1 -HsT -sT)/(2*sqrt(s*s1)) ) 
      -  Hz * log( K/(2.*sqrt(s*s1)*(sT-z) ) )
      ) / (16.*M_PI*M_PI*s) , 0 );
}

Expression AmpGen::phaseSpace(const Expression& s, const Particle& p, const size_t& l)
{
  auto fs = p.getFinalStateParticles();
  if( fs.size() == 2 ){ 
    auto s1 = p.daughter(0)->massSq();
    auto s2 = p.daughter(1)->massSq();
    auto k2 = make_cse( Q2(s,s1,s2) );
    const Expression k2p  = Ternary( k2 > 0, k2, 0 );
    auto phsp_parameterisation = p.attribute("phsp"); 
    if( !phsp_parameterisation ){
      const Expression radius       = Parameter(p.name()  + "_radius", p.props()->radius());
      return rho_twoBody(s, s1, s2) * BlattWeisskopf(k2p*radius*radius, l);
    }
    if( phsp_parameterisation == std::string("arXiv.0707.3596") ){
      INFO("Got AS parametrisation");
      Expression E = 1;
 //     return 2*fpow(k2,l)*complex_sqrt(k2) * E * 2 / ((1.5+k2p) * sqrt(s) ); 
      return 2 * complex_sqrt(k2/s);
    }
  }
  if( fs.size() == 3 && ! p.daughter(0)->isStable() ) return rho_threeBody( s, *p.daughter(0), *p.daughter(1) );
  if( fs.size() == 3 && ! p.daughter(1)->isStable() ) return rho_threeBody( s, *p.daughter(1), *p.daughter(0) );
  ERROR("Phase space only implemented for two and three body decays");
  return 0;
}

DEFINE_LINESHAPE( CoupledChannel )
{
  const auto props              = ParticlePropertiesList::get( particleName );
  const Expression mass         = Parameter( particleName + "_mass" , props->mass() );
  const Expression width        = Parameter( particleName + "_width", props->width() );
  const Expression radius       = Parameter( particleName + "_radius", props->radius() );
  const Expression I = Constant(0.,1.);
  std::vector<std::string> channels = NamedParameter<std::string>( particleName + "_channels").getVector();  
  Expression totalWidth = 0; 
  Expression totalWidthAtPole = 0 ; 
  ADD_DEBUG( s , dbexpressions );
  for( size_t i = 0 ; i < channels.size(); i+=2 ){
    Particle p( channels[i] ); 
    DEBUG( "Adding channel ... " << p.uniqueString() << " coupling = " << NamedParameter<std::string>( channels[i+1]  ) );
    Expression coupling = Parameter(channels[i+1], 0);
    totalWidth       = totalWidth       + coupling * phaseSpace(s        , p, p.orbital());
    totalWidthAtPole = totalWidthAtPole + coupling * phaseSpace(mass*mass, p, p.orbital());    
    ADD_DEBUG(coupling, dbexpressions);
    ADD_DEBUG(phaseSpace(s,p,p.orbital() ), dbexpressions);
    ADD_DEBUG(phaseSpace(mass*mass,p,p.orbital() ), dbexpressions);
  }
  ADD_DEBUG(totalWidth, dbexpressions);
  ADD_DEBUG(totalWidthAtPole, dbexpressions);
  const Expression q2  = make_cse(abs(Q2(s, s1, s2)));
  Expression formFactor                        = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( lineshapeModifier == "BL" )  formFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  const Expression widthNorm = lineshapeModifier == "norm" ? width / totalWidthAtPole : 1; 
  const Expression D  = mass*mass - Constant(0,1)*mass*totalWidth*widthNorm; 
  const Expression BW = 1. / ( D - s ); 
  const Expression kf = kFactor( mass, width, dbexpressions );
  ADD_DEBUG( totalWidth , dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return kf * formFactor * BW;
}
