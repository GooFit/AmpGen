#ifndef AMPGEN_DYNAMICFCN_H
#define AMPGEN_DYNAMICFCN_H

#include "AmpGen/MsgService.h"
#include <dlfcn.h>
#include <iostream>

namespace AmpGen
{
  template<typename RETURN_TYPE>
  class DynamicFCN;

  template <typename RETURN_TYPE, typename... IN_TYPES>
  class DynamicFCN<RETURN_TYPE( IN_TYPES... )>
  {

  private:
    void* m_handle;
    RETURN_TYPE ( *m_fcn )( IN_TYPES... );

  public:
    DynamicFCN() : m_handle( nullptr ), m_fcn( nullptr ) {}
    DynamicFCN( const std::string& lib, const std::string& name ) :

      m_handle(dlopen( lib.c_str(), RTLD_NOW )) {
        set(m_handle,name);
    }
    DynamicFCN( void* handle, const std::string& name ) : m_handle(handle) { set( handle, name ); }
    
    bool set( const std::string& lib, const std::string& name )
    {
      DEBUG("Linking handle: " << lib << ":" << name );
      m_handle = dlopen( lib.c_str(), RTLD_NOW );
      if( !m_handle ){
        DEBUG( dlerror() );
        return false;
      }
      return set(m_handle,name);
    }
    bool set( void* handle, const std::string& name )
    {
      m_fcn = (RETURN_TYPE( * )( IN_TYPES... ))dlsym( handle, name.c_str() );
      if ( m_fcn == nullptr ) {
        ERROR( dlerror() );
        return false;
      }
      return true;
    }
    RETURN_TYPE operator()( IN_TYPES... input ) const { return ( *m_fcn )( input... ); }
    bool isLinked() const { return m_fcn != nullptr; }
    ~DynamicFCN() = default ;
  };
} // namespace AmpGen

#endif
