#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

DEFINE_LINESHAPE( Poly )
{

  size_t degree = NamedParameter<size_t>( lineshapeModifier + "::Degree" );
  auto params   = parameterVector(lineshapeModifier + "_c", degree + 1);
  ADD_DEBUG(s, dbexpressions);
  if( dbexpressions != nullptr ) for( size_t i=0;i<params.size();++i ) 
    dbexpressions->emplace_back( lineshapeModifier +"_c"+std::to_string(i), params[i] );

  return pol(s, params);
}
