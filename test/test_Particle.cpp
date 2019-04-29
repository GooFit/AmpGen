#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Particle"

#include <boost/test/unit_test.hpp>
#include "AmpGen/Particle.h"

BOOST_AUTO_TEST_CASE( Particle_quasiStableTree )
{
  auto d1 = AmpGen::Particle("D~0{D0{K*(892)0{K+,pi-},rho(770)0{pi+,pi-}}}");
  auto d2 = AmpGen::Particle("B+{D0{K*(892)0{K+,pi-},rho(770)0{pi+,pi-}},K+}");
  auto d3 = AmpGen::Particle("B+{D~0{D0{K*(892)0{K+,pi-},rho(770)0{pi+,pi-}}},K+}");
  BOOST_CHECK( d1.quasiStableTree().decayDescriptor() == "D~0{K+,pi+,pi-,pi-}" );
  BOOST_CHECK( d2.quasiStableTree().decayDescriptor() == "B+{D0{K+,pi+,pi-,pi-},K+}" );
  BOOST_CHECK( d3.quasiStableTree().decayDescriptor() == "B+{D~0{K+,pi+,pi-,pi-},K+}" );
}


