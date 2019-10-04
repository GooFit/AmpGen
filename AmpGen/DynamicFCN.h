#ifndef AMPGEN_DYNAMICFCN_H
#define AMPGEN_DYNAMICFCN_H

#include "AmpGen/MsgService.h"
#include <dlfcn.h>
#include <iostream>

namespace AmpGen
{
  /**@class DynamicFCN 
     @brief Wrapper to give templated interface to a function contained in a dynamically linked library. 
    @tparam RETURN_TYPE type that this function returns
    @tparam IN_TYPES input types for this function 
    DynamicFCN is a template wrapper for a C-style function pointer and a handle to a dynamic library, i.e. provides a more C++-like interface to dlsym. 
   For example, for a function called "foo" in library "bar.so" that returns a double from two doubles can be used as
    
     \code{cpp}
       DynamicFCN<double(double,double)> foo("bar.so","foo");
       double z = foo(4,5); 
     \endcode  
  
     i.e. the foo function object can then be used as a regular function.*/

  template <class RETURN_TYPE, class ...IN_TYPES> 
  class DynamicFCN; 

  template <class RETURN_TYPE, class ...IN_TYPES>
  class DynamicFCN<RETURN_TYPE( IN_TYPES... )>
  {
  private:
    void* m_handle                        = {nullptr};
    RETURN_TYPE ( *m_fcn )( IN_TYPES... ) = {nullptr};

  public:
    DynamicFCN() = default; 
    DynamicFCN( const std::string& lib, const std::string& name ) : 
      m_handle(dlopen( lib.c_str(), RTLD_NOW )) {
        set(m_handle,name);
    }
    DynamicFCN( void* handle, const std::string& name ) : m_handle(handle) { set( handle, name ); }
    ~DynamicFCN() = default;

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
  };
} // namespace AmpGen

#endif
