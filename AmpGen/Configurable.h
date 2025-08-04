#ifndef AMPGEN_CONFIGURABLE_H 
#define AMPGEN_CONFIGURABLE_H 1 
#include "AmpGen/MsgService.h" 

namespace AmpGen { 
  class Configurable { 
    public: 
      void registerParameter( const std::string& name ){
        INFO("Registering: "<< name ); 
        m_parameters.push_back( name ); 
      }
    private: 
      std::vector<std::string> m_parameters; 
  }; 
}
#endif 
