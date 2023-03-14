#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Particle.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/CoupledChannel.h"
#include "AmpGen/Pade.h"
#include "AmpGen/NumericalIntegration.h"

using namespace AmpGen;
using namespace AmpGen::fcn; 
using namespace std::complex_literals; 

// ENABLE_DEBUG( Lineshape::CoupledChannel );

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

namespace {
  template <typename fcn_type> auto analytic_continuation( fcn_type&& functor, const double& s_thresh )
  { 
    return [functor, s_thresh]( const double& s )
    {
      double epsilon=1e-9;
      auto real = (s-s_thresh) * integrate_1d_cauchy( [=](const double& s){ return functor(s) / (s-s_thresh) ; } , s, s_thresh +  epsilon, 100000 ) / M_PI; 
      return 16 * M_PI * std::complex<double>( real, s > s_thresh ? functor(s) : 0 ); 
    };
  }

  template <typename fcn_type> auto real( fcn_type&& functor )
  {
    return [functor]( const double& s){ return std::real(functor(s)); }; 
  }

  template <typename fcn_type> auto imag( fcn_type&& functor )
  {
    return [functor]( const double& s){ return std::imag(functor(s)); }; 
  }

  template <unsigned L, typename expression_type> expression_type blatt_weisskopf_sq( const expression_type& z )
  {
    if constexpr( L == 0 )      return expression_type(1);
    else if constexpr( L == 1 ) return expression_type( 2 * z / ( 1 + z )  ); 
    else if constexpr( L == 2 ) return expression_type( 13 * z * z / ( z*z + 3*z * 9 ) ); 
    else if constexpr( L == 3 ) return expression_type( 277  * z * z * z  / (z*(z-15)*(z-15) + 9* ( 2 * z - 5 ) * (2*z-5 ) ) );  
    else {
      std::cout << "error, blatt weisskopf not implemented for L> 3 (if you are trying to calculate dispersive corrections in such a case, good luck to ya" << std::endl; 
      return 0; 
    }
  }
  template <typename T, typename T2>
  T q2( const T& s, const T2& s1, const T2& s2 )
  {
    return 0.25 * (s - 2*(s1+s2) + ( s1-s2)*(s1-s2)/s );
  }

  template <unsigned N> Expression rho_pade(const Expression& s, const Particle& p )
  {
    double m1 = p.daughter(0)->props()->mass();
    double m2 = p.daughter(1)->props()->mass();
    double s_thresh = (m1+m2)*(m1+m2);
    double R  = p.props()->radius(); 
    auto CM_num = analytic_continuation( [=](const double& s){ 
      return blatt_weisskopf_sq<N>(q2(s,m1*m1,m2*m2) * R *R ) * sqrt( q2(s,m1*m1,m2*m2) / s ) / ( 8 * M_PI ); }, s_thresh ); 
    Pade <4, double> cm_pade_lo( real(CM_num), 0.0, s_thresh );
    Pade <4, double> cm_pade_hi( real(CM_num), s_thresh, 5 );
    return Ternary( s >= s_thresh, cm_pade_hi(s) + 1i * blatt_weisskopf_sq<N>(q2(s,m1*m1,m2*m2) * R * R ) * 2 * sqrt( Q2(s,m1*m1,m2*m2) / s ),  cm_pade_lo(s) )/ 1i;
  }
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
    else if( phsp_parameterisation == std::string("arXiv.0707.3596") ){
      INFO("Got AS parametrisation");
      return 2 * complex_sqrt(k2/s);
    }
    else if( phsp_parameterisation == std::string("CM") )
    {
      if( l == 1 ) return rho_pade<1>(s, p); 
      if( l == 2 ) return rho_pade<2>(s, p); 
      if( l == 3 ) return rho_pade<3>(s, p); 
      auto m1 = p.daughter(0)->mass();
      auto m2 = p.daughter(1)->mass();
      Expression s1 = m1*m1;
      Expression s2 = m2*m2;
      Expression sT = (m1+m2)*(m1+m2);
      Expression q2 = (s*s + s1*s1 + s2*s2 - 2*s*s1 - 2*s*s2 - 2*s1*s2)/(4*s);
      auto q   = fcn::complex_sqrt(q2);
      auto arg = (s1 + s2 - s + 2*fcn::sqrt(s) * q )/ (2*m1*m2);
      return (2.*q * fcn::log(arg) / fcn::sqrt(s) - (s1-s2)*( 1./s - 1./sT ) * fcn::log(m1/m2) ) / ( 16.i * M_PI * M_PI ); 
    }
    else {
      FATAL("Parametrisation: " << *phsp_parameterisation << " not found");
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
    totalWidth       += coupling * phaseSpace(s        , p, p.L());
    totalWidthAtPole += coupling * phaseSpace(mass*mass, p, p.L());    
    ADD_DEBUG(coupling, dbexpressions);
    ADD_DEBUG(phaseSpace(s,p,p.L() ), dbexpressions);
    ADD_DEBUG(phaseSpace(mass*mass,p,p.L() ), dbexpressions);
  }
  ADD_DEBUG(totalWidth, dbexpressions);
  ADD_DEBUG(totalWidthAtPole, dbexpressions);
  const Expression q2  = make_cse(fcn::abs(Q2(s, s1, s2)));
  Expression formFactor                        = fcn::sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if ( lineshapeModifier == "BL" )  formFactor = fcn::sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  const Expression widthNorm = lineshapeModifier == "norm" ? width / totalWidthAtPole : 1; 
  const Expression D  = mass*mass - Constant(0,1)*mass*totalWidth*widthNorm; 
  const Expression BW = 1. / ( D - s ); 
  const Expression kf = kFactor( mass, width, dbexpressions );
  ADD_DEBUG( totalWidth , dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return kf * formFactor * BW;
}
