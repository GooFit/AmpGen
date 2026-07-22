#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticlePropertiesList.h"

using namespace AmpGen;

DEFINE_LINESHAPE(Photon) {
  auto props = ParticlePropertiesList::get(particleName);
  Expression mass = Variable(particleName + "_mass", props->mass());
  return 1. / (s - mass * mass);
}
