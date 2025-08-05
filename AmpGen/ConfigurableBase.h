#ifndef AMPGEN_CONFIGURABLEBASE_H
#define AMPGEN_CONFIGURABLEBASE_H 1
#include "AmpGen/MsgService.h"
#include <vector> 
#include <string> 

namespace AmpGen {
  class ConfigurableBase { 
    public:
    void registerParameter(const std::string &name) {
      m_parameters.push_back(name);
    }

    void print() const {
      for( auto const& p : m_parameters ) INFO( p ); 
    }

    std::string remove_namespace( const std::string& name ); 
  private:
    std::vector<std::string> m_parameters;
  };
}
#endif
