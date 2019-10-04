#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "AmplitudeRules"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"

using namespace AmpGen; 

BOOST_AUTO_TEST_CASE( test_AmplitudeRule ) {

  MinuitParameter re = MinuitParameter("D0{K*(892)bar0{K-,pi+},pi0}_Re", Flag::Free,1.,0.);
  MinuitParameter im = MinuitParameter("D0{K*(892)bar0{K-,pi+},pi0}_Im", Flag::Free,0.,0.);

  AmplitudeRule test(&re,&im);

  BOOST_CHECK( test.name()   == "D0{K*(892)bar0{K-,pi+},pi0}" );
  BOOST_CHECK( test.head()   == "D0");
  BOOST_CHECK( test.prefix() == "" );
  BOOST_CHECK( test.eventType() == EventType({"D0","K-","pi+","pi0"}) );
}
AmplitudeRules rule_set( 
    MinuitParameterSet( {
   new MinuitParameter("D0{K*(892)bar0,pi0}_Re", Flag::Fix,1.,0.)
 , new MinuitParameter("D0{K*(892)bar0,pi0}_Im", Flag::Fix,1.,0.)
 , new MinuitParameter("D0{rho(770)+,pi-}_Re"  , Flag::Fix,1.,0.)
 , new MinuitParameter("D0{rho(770)+,pi-}_Im"  , Flag::Fix,2.,0.)
 , new MinuitParameter("K*(892)bar0{K-,pi+}_Re", Flag::Free,sqrt(1./3.),0.)
 , new MinuitParameter("K*(892)bar0{K-,pi+}_Im", Flag::Free,0,0.)
 , new MinuitParameter("K*(892)bar0{K0,pi0}_Re", Flag::Fix,sqrt(2./3.),0.)
 , new MinuitParameter("K*(892)bar0{K0,pi0}_Im", Flag::Fix,0,0.) } ) ); 

BOOST_AUTO_TEST_CASE( test_AmplitudeRules_constructor ){
  BOOST_CHECK( rule_set.rules().size() == 2 ); /// number of head decays
  BOOST_CHECK( rule_set.hasDecay("D0") == true ); /// has decays for D0 
  BOOST_CHECK( rule_set.hasDecay("a(1)(1260)+") == false ); /// has decays for D0 
}

BOOST_AUTO_TEST_CASE( test_couplingConstant )
{
  auto matches = rule_set.getMatchingRules( EventType({"D0","K-","pi+","pi0"}) );  
  auto matches2 = rule_set.getMatchingRules( EventType({"D0","K0","pi0","pi0"}) ); 
  
  BOOST_TEST( std::real(matches[0].second()) == 1./sqrt(3.), boost::test_tools::tolerance(1e-10) );  
  BOOST_TEST( std::imag(matches[0].second()) == 1./sqrt(3.), boost::test_tools::tolerance(1e-10) );  
  BOOST_TEST( matches[0].second.isFixed() == false );
  BOOST_TEST( matches2[0].second.isFixed() == true );
  BOOST_TEST( matches[0].second.contains("K*(892)bar0{K0,pi0}") == false );
  BOOST_TEST( matches2[0].second.contains("K*(892)bar0{K0,pi0}") == true );
}

