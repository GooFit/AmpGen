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
  std::vector<std::string> forbidden = NamedParameter<std::string>("AddCPConjugate::Forbid").getVector();  
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
        auto test_particle = Particle(pname);
        if( std::count( forbidden.begin(), forbidden.end(), test_particle.name() ) != 0 ){
          INFO( "Name:" << test_particle.name() << " is forbidden" );
          continue; 
        }  
        Particle test = Particle(test_particle).conj();
        if( cartOrPolar == "polar" )     sgn = reOrIm == "Re" ? test.CP() : 1; 
        if( cartOrPolar == "cartesian" ) sgn = test.CP();
        new_name = test.uniqueString() +"_"+reOrIm;
      }
      else if( tokens.size() == 2 ) {
        auto props = AmpGen::ParticlePropertiesList::get(pname, true);
        if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
      }
    }
    if( mps.find( new_name ) == nullptr )
    {
      INFO("Adding CP conjugated rule: " << new_name << " from " << name ); 
      tmp.push_back( new MinuitExpression(new_name, sgn * MinuitParameterLink(param) )) ;  
    }
  }
  for( auto& p : tmp ) mps.add( p );
}
