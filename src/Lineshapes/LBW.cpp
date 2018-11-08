#include <stddef.h>
#include <string>
#include <vector>

#include "AmpGen/Lineshapes.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/ParticleProperties.h"

using namespace AmpGen;
using namespace AmpGen::fcn;

Expression width_no_norm( const Expression& s, 
                          const Expression& q2, 
                          const Expression& radius,    
                          const unsigned int& L)
{
  Expression BF        = BlattWeisskopf_Norm( q2 * radius * radius, 0, L );
  return BF * sqrt( q2 ) * fpow( q2,L ) / sqrt(s);
}

DEFINE_LINESHAPE( LBW )
{
  auto props = ParticlePropertiesList::get( particleName );
  const Expression& mass     = Parameter( particleName + "_mass", props->mass() );
  const Expression& width0   = Parameter( particleName + "_width", props->width() );
  const Expression& radius   = Parameter( particleName + "_radius", props->radius() );  
  auto  orbital_parameters   = NamedParameter<size_t>( particleName +"_waves").getVector();
  const Expression s0        = mass*mass;  
  const Expression q2_sgned  = make_cse( Q2( s , s1, s2 ) ) ;
  const Expression q20       = make_cse( Q2( s0, s1, s2 ) );
  const Expression J = Constant(0,1);
  const Expression q2  = Ternary( q2_sgned > 0, q2_sgned, 0 );

  Expression width_total = 0;
  Expression width_norm  = 0; 
  for( auto& o : orbital_parameters ){
    auto g = Parameter( particleName + "_g"+std::to_string(o) );
    width_total = width_total + g * width_no_norm(s ,q2 ,radius, o);
    width_norm  = width_norm  + g * width_no_norm(s0,q20,radius, o);
    ADD_DEBUG( g , dbexpressions );
    ADD_DEBUG( width_total, dbexpressions );
  }
  const Expression FormFactor   = fcn::sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) ); 
  const Expression runningWidth = width0 * width_total / width_norm ; 
  const Expression BW           = FormFactor / ( s0 - s  -J*mass * runningWidth );
  const Expression kf           = kFactor( mass, width0, dbexpressions );
  ADD_DEBUG( FormFactor, dbexpressions );
  ADD_DEBUG( runningWidth, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  ADD_DEBUG( kf, dbexpressions );
  return kf * BW;
}
