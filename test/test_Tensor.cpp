#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Tensor"

#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "AmpGen/Tensor.h"

using namespace AmpGen; 

BOOST_AUTO_TEST_CASE ( constructor )
{
  Tensor A({"A_11","A_12","A_13",
      "A_21","A_22","A_23"}, Tensor::dim(2,3) );


}

BOOST_AUTO_TEST_CASE( multiply )
{
  Tensor::Index a,b,c;
  Tensor A({1,2,3,
      4,5,6},Tensor::dim(2,3) );
  Tensor B({7,8,
      9,10,
      11,12}, Tensor::dim(3,2) );

  Tensor C = A(a,b)*B(b,c);

  BOOST_TEST( int( std::real( C[{0,0}]()) ) == 58 );
  BOOST_TEST( int( std::real( C[{0,1}]()) ) == 64 );
  BOOST_TEST( int( std::real( C[{1,0}]()) ) == 139 );
  BOOST_TEST( int( std::real( C[{1,1}]()) ) == 154 );

}

BOOST_AUTO_TEST_CASE( inverse )
{
  Tensor t( Tensor::dim(5,5) );
  t(0, 0) = 0.999741748906672;
  t(0, 1) = 0.162909875391051;
  t(0, 2) = 0.28261780529283;
  t(0, 3) = 0.947201082017273;
  t(0, 4) = 0.23165654274635;
  t(1, 0) = 0.484973614336923;
  t(1, 1) = 0.957476956536993;
  t(1, 2) = 0.744305343134329;
  t(1, 3) = 0.540043658344075;
  t(1, 4) = 0.739952981472015;
  t(2, 0) = 0.759943798184395;
  t(2, 1) = 0.65863661444746;
  t(2, 2) = 0.315637621562928;
  t(2, 3) = 0.804403014713898;
  t(2, 4) = 0.519672115100548;
  t(3, 0) = 0.168572421884164;
  t(3, 1) = 0.475529730319977;
  t(3, 2) = 0.392313994001597;
  t(3, 3) = 0.221667687175795;
  t(3, 4) = 0.213190458947793;
  t(4, 0) = 0.0303352042101324;
  t(4, 1) = 0.333539249841124;
  t(4, 2) = 0.194148850627244;
  t(4, 3) = 0.943716780748218;
  t(4, 4) = 0.579931674990803;

  auto it = t.Invert();

  BOOST_TEST( std::real( it(0, 0)() ) == 0.533179172025328, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 1)() ) == 0.754080023420741, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 2)() ) == 0.564168881565736, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 3)() ) == -1.75645898717539, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 4)() ) == -1.03498527412253, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 0)() ) == -1.75316998902971, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 1)() ) == -2.33293061410697, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 2)() ) == 3.05756349746262, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 3)() ) == 3.38025891202066, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 4)() ) == -0.305513880783314, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 0)() ) == 1.77658898428149, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 1)() ) == 1.51983750058763, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 2)() ) == -3.57990876398275, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 3)() ) == 1.20923813570621, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 4)() ) == 0.114520342807423, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 0)() ) == 0.282573916596254, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 1)() ) == -1.73119940397175, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 2)() ) == 0.151145803840974, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 3)() ) == 2.43137060997662, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 4)() ) == 1.06677185399228, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 0)() ) == -0.0741734831238967, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 1)() ) == 3.61065979016023, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 2)() ) == -0.835504067045532, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 3)() ) == -6.21360079842779, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 4)() ) == 0.179905447460274, boost::test_tools::tolerance(1e-10) );

}


BOOST_AUTO_TEST_CASE( inverse_symmetric )
{
  Tensor t_symmetric( Tensor::dim(5,5) );
  t_symmetric(0, 0) = 0.898304857779294;
  t_symmetric(0, 1) = t_symmetric(1, 0) = 0.665563931455836;
  t_symmetric(0, 2) = t_symmetric(2, 0) = 0.49861030979082;
  t_symmetric(0, 3) = t_symmetric(3, 0) = 0.560628257226199;
  t_symmetric(0, 4) = t_symmetric(4, 0) = 0.182284645503387;
  t_symmetric(1, 1) = 0.296525530749932;
  t_symmetric(1, 2) = t_symmetric(2, 1) = 0.117408933350816;
  t_symmetric(1, 3) = t_symmetric(3, 1) = 0.0629176658112556;
  t_symmetric(1, 4) = t_symmetric(4, 1) = 0.648125574691221;
  t_symmetric(2, 2) = 0.725418528541923;
  t_symmetric(2, 3) = t_symmetric(3, 2) = 0.637131158262491;
  t_symmetric(2, 4) = t_symmetric(4, 2) = 0.713885062374175;
  t_symmetric(3, 3) = 0.0995762431994081;
  t_symmetric(3, 4) = t_symmetric(4, 3) = 0.699267196236178;
  t_symmetric(4, 4) = 0.107812469825149;

//  t_symmetric.imposeSymmetry(0,1);
  auto it = t_symmetric.Invert();

  BOOST_TEST( std::real( it(0, 0)() ) == 1.66189984270993, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 1)() ) == -0.688951871462912, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 0)() ) == -0.688951871462912, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 2)() ) == 0.148581580245623, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 0)() ) == 0.148581580245622, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 3)() ) == 0.276468177065467, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 0)() ) == 0.276468177065467, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(0, 4)() ) == -1.44516486843412, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 0)() ) == -1.44516486843412, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 1)() ) == 1.60872622199432, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 2)() ) == -1.51350082633105, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 1)() ) == -1.51350082633105, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 3)() ) == -0.04273326919781, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 1)() ) == -0.0427332691978104, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(1, 4)() ) == 1.79270863535293, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 1)() ) == 1.79270863535293, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 2)() ) == -0.850132539931173, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 3)() ) == 2.15728970691434, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 2)() ) == 2.15728970691434, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(2, 4)() ) == 0.484447195530136, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 2)() ) == 0.484447195530137, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 3)() ) == -2.16626046580015, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(3, 4)() ) == -0.444859428711916, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 3)() ) == -0.444859428711916, boost::test_tools::tolerance(1e-10) );
  BOOST_TEST( std::real( it(4, 4)() ) == 0.619288662427455, boost::test_tools::tolerance(1e-10) );

}
