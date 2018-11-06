#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "PhaseSpace"

#include <boost/test/unit_test.hpp>
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/CompiledExpression.h"
#include "TRandom3.h"

namespace tt = boost::test_tools;
using namespace AmpGen;

BOOST_AUTO_TEST_CASE( phaseSpace_threeBody )
{
  AmpGen::PhaseSpace phsp( AmpGen::EventType( {"Xi(c)+","Lambda0","K0","pi+"} ) , new TRandom3(15) );

  std::vector<double> test_event = {-0.235918, -0.242689, 0.278177, 1.19862,-0.300608, -0.0584944, -0.0117436, 0.584418, 0.536526, 0.301183, -0.266434, 0.684864} ;
  
  auto new_event = phsp.makeEvent(0);

  for( int i = 0 ; i < 12 ; ++i )
    BOOST_TEST( test_event[i] == new_event[i]  , boost::test_tools::tolerance(1e-5) );
}

