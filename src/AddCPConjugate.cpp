#include <string>
#include <vector>

#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Property.h"
#include "AmpGen/ExpressionParser.h"
#include "AmpGen/MinuitExpression.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AmplitudeRules.h" 

using namespace AmpGen; 


void AmpGen::AddCPConjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  auto coord = Coupling(nullptr, nullptr).coordinates(); 
  std::vector<std::string> forbidden = Property<std::vector<std::string>>(nullptr, "AddCPConjugate::Forbid"); 

  for( auto& param : mps ){
    const std::string name = param->name();
    if( name.find("::") != std::string::npos ){
      auto pos = name.find("::");
      auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
      if( props != 0 )
      {
        std::string new_name = props->anti().name() + name.substr(pos); 
        tmp.push_back( new MinuitExpression(new_name, MinuitParameterLink(param) )) ;  
      }
      continue;
    }
    auto tokens=split(name,'_');
    
    if( tokens.size() == 1 && Particle::isValidDecayDescriptor(tokens[0]) )
    {
      auto particle = Particle(tokens[0]); 
      if( std::count( forbidden.begin(), forbidden.end(), particle.name() ) == 0 ){
        //DEBUG( particle.name() << "->" << particle );   [ to do: fix DEBUG for free functions ] 
        tmp.push_back( new MinuitExpression( particle.conj().decayDescriptor(), fcn::conj( MinuitParameterLink(param)) ) );
      }
      else {
        //DEBUG("Skipping: " << particle ); 
      }
      continue; 
    }
    std::string reOrIm = *tokens.rbegin();
    std::string pname   = tokens[0];
    std::string new_name = ""; 
    int sgn=1;
    if ( reOrIm == "Re" || reOrIm == "Im" ){
      auto test_particle = Particle(pname);
      //DEBUG( test_particle.name() << "->" << test_particle ); 
      if( std::count( forbidden.begin(), forbidden.end(), test_particle.name() ) != 0 ){
        //DEBUG("Skipping: " << test_particle ); 
        continue; 
      } 
      Particle test = Particle(test_particle).conj();
      if( ! test.props()->hasDistinctAnti() ) continue; 
      if( coord == coordinateType::polar )     sgn = reOrIm == "Re" ? test.CP() : 1; 
      if( coord == coordinateType::cartesian ) sgn = test.CP();
      new_name = test.uniqueString() +"_"+reOrIm;
    }
    else if( tokens.size() == 2 ) {
      auto props = AmpGen::ParticlePropertiesList::get(pname, true);
      if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
    }
    tmp.push_back( new MinuitExpression(new_name, sgn * MinuitParameterLink(param) )) ;  
  }
  for( auto& p : tmp ){
    if( mps.find(p->name()) == 0 ) mps.add( p );
    else delete p; 
  }
}
