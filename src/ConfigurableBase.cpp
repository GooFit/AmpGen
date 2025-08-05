#include "AmpGen/ConfigurableBase.h" 
#include "AmpGen/Utilities.h" 

using namespace AmpGen; 

std::string ConfigurableBase::remove_namespace ( const std::string& name ){
  return replaceAll( name, "AmpGen::", ""); 
} 
