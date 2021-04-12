#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

DEFINE_GENERIC_SHAPE( PolyNR )
{
  auto p0 = *p.daughter(0);
  auto p1 = *p.daughter(1);
  auto p2 = *p.daughter(2);
  auto p01 = (p0.P() + p1.P());
  auto p02 = (p0.P() + p2.P());
  auto s01 = dot( p01, p01 );
  auto s02 = dot( p02, p02 );
  size_t degree = NamedParameter<size_t>( lineshapeModifier + "::Degree" ) +1;
  std::vector< std::vector<Parameter> > C( degree, std::vector<Parameter>(degree) );
  
  for( size_t i = 0 ; i != degree ; ++i )
  {
    for( size_t j = 0 ; j != degree ; ++j ){
      auto pname = lineshapeModifier +"_"+std::to_string(i) + "_"+std::to_string(j) ;
      C[i][j] = Parameter(pname, 0);
      if( dbexpressions != nullptr )
        dbexpressions->emplace_back( pname, C[i][j] );
    }
  }
  Expression rt; 
  for( size_t i = 0 ; i != degree; ++i )
  {
    for( size_t j = 0 ; j != degree ; ++j )
    {
      rt = rt + C[i][j] * fcn::pow(s01, i) * fcn::pow(s02, j); 
    }
  }
  return rt;
}
