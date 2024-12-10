#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Expression"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"

#if INSTRUCTION_SET != 0 
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

using namespace AmpGen;
BOOST_AUTO_TEST_CASE ( test_constructor )
{
  MinuitParameterSet mps;
  mps.add( "D0{K0S0,rho(770)0{pi+,pi-}}_Re", Flag::Fix, 1, 0); 
  mps.add( "D0{K0S0,rho(770)0{pi+,pi-}}_Im", Flag::Fix, 0, 0); 
  mps.add( "D0{K*(892)+{K0S0,pi+},pi-}_Re", Flag::Fix, 0.5, 0); 
  mps.add( "D0{K*(892)+{K0S0,pi+},pi-}_Im", Flag::Fix, 0.2, 0); 
  auto t = EventType({"D0","K0S0","pi+","pi-"});
  
  auto events = EventList_type(Generator<>(t).generate(500)); 
  
  CoherentSum test( t, mps );
  test.prepare();
  auto eval = test.amplitudeEvaluator( &events ); 

  for( auto& event : events )
  { 
    BOOST_TEST( eval(event) == test.getValNoCache(event), boost::test_tools::tolerance(1e-14) ); 
  }
}
