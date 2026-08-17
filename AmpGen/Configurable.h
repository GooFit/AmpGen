#ifndef AMPGEN_CONFIGURABLE_H
#define AMPGEN_CONFIGURABLE_H 1
#include "AmpGen/MsgService.h"
#include "AmpGen/Property.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/ConfigurableBase.h"

namespace AmpGen {

  template <typename T> class Configurable : public ConfigurableBase {
    static T *gImpl;

  public:
    virtual ~Configurable() = default;

  protected:
    Property<bool> m_verbose{this, remove_namespace(type_string<T>()) + "::Verbose", false, "Enable verbose printing"};
  };
}
#endif
