#include <string>
#include <vector>

#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ExpressionParser.h"
#include "AmpGen/MinuitExpression.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticlePropertiesList.h"

using namespace AmpGen; 

void AmpGen::AddCPConjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  std::string cartOrPolar = NamedParameter<std::string>("CouplingConstant::Coordinates" ,"cartesian");
  for( auto& param : mps ){
    const std::string name = param->name();
    std::string new_name = name; 
    int sgn=1;
    if( name.find("::") != std::string::npos ){
      auto pos = name.find("::");
      auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
      if( props != 0 ) new_name = props->anti().name() + name.substr(pos); 
    }
    else { 
      auto tokens=split(name,'_');
      std::string reOrIm = *tokens.rbegin();
      std::string pname   = tokens[0];
      if ( reOrIm == "Re" || reOrIm == "Im" ){
        Particle test = Particle(pname).conj();
        if( cartOrPolar == "polar" )     sgn = reOrIm == "Re" ? test.CP() : 1; 
        if( cartOrPolar == "cartesian" ) sgn = test.CP();
        new_name = test.uniqueString() +"_"+reOrIm;
      }
      else if( tokens.size() == 2 ) {
        auto props = AmpGen::ParticlePropertiesList::get(name, true);
        if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
      }
    }
    if( mps.find( new_name ) == nullptr )
    {
      tmp.push_back( new MinuitExpression(new_name, sgn * MinuitParameterLink(param) )) ;  
    }
  }
  for( auto& p : tmp ) mps.add( p );
}
