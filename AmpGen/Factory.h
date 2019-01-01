#ifndef AMPGEN_FACTORY_H
#define AMPGEN_FACTORY_H

#include <cxxabi.h>
#include <map>

#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

#define REGISTER( BASE_CLASS, DERIVED_CLASS )                                                                          \
  std::string DERIVED_CLASS::_id = AmpGen::Factory<BASE_CLASS>::Register( #DERIVED_CLASS, new DERIVED_CLASS() )
#define REGISTER_WITH_KEY( BASE_CLASS, DERIVED_CLASS, KEY, KEY_TYPE )                                                  \
  KEY_TYPE DERIVED_CLASS::_id = AmpGen::Factory<BASE_CLASS, KEY_TYPE>::Register( KEY, new DERIVED_CLASS() )

namespace AmpGen
{
  /**@class Factory
   * Static factory to construct classes from a hierarchy based on a key (normally std::string)
   */
  template <class TYPE, class KEY_TYPE = std::string>
  class Factory
  {
  public:
    std::map<KEY_TYPE, TYPE*> m_terms;     ///< map of objected made by this factory
    static Factory<TYPE, KEY_TYPE>* gImpl; ///< pointer to static implementation

    static Factory<TYPE, KEY_TYPE>* getMe()
    {
      if ( !gImpl ) gImpl = new Factory<TYPE, KEY_TYPE>();
      return gImpl;
    }
    static TYPE* get( const KEY_TYPE& type, const bool quiet = false )
    {
      auto ptrToStatic = getMe();
      auto raw_base    = ptrToStatic->m_terms.find( type );
      if ( raw_base == ptrToStatic->m_terms.end() ) {
        if ( !quiet ) ERROR( type << " not found in Factory<" << typeof<TYPE>() << typeof<KEY_TYPE>() << " >" );
        return nullptr;
      }
      auto objectToReturn = raw_base->second->create();
      return objectToReturn;
    }
    static KEY_TYPE Register( const KEY_TYPE& key, TYPE* object )
    {
      getMe()->m_terms[key] = object;
      return key;
    }
  };

  template <class TYPE, class KEY_TYPE> 
    Factory<TYPE,KEY_TYPE>* Factory<TYPE,KEY_TYPE>::gImpl = nullptr; 
} // namespace AmpGen

#endif

