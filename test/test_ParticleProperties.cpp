#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "ParticlePropertiesList"

#include <boost/test/unit_test.hpp>
#include "AmpGen/ParticlePropertiesList.h"

BOOST_AUTO_TEST_CASE( PDG_mass )
{
  BOOST_CHECK ( AmpGen::ParticlePropertiesList::get("K*(892)0")->mass() == 0.89600 );
}


