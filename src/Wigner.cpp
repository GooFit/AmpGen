#include "AmpGen/Wigner.h"

#include <memory.h>
#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <complex>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Simplify.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Vertex.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;
using namespace AmpGen::fcn;
using namespace std::complex_literals; 

double fact( const double& z )
{
  double f=1;
  for( int i=1;i<=z;++i) f*=i;
  return f;
}

double nCr_( const int& n, const int& r ){
  double z=1;
  for( int f=1; f <= r ; ++f ) z *= double(n+1-f)/double(f);
  return z;
}

Expression expandedBinomial( const Expression& x, const unsigned int& n )
{
  Expression sum;
  for( unsigned int k = 0 ; k <= n ; ++k ) sum = sum + nCr_(n,k) * fcn::fpow(x,k);
  return sum; 
}

Expression AmpGen::wigner_d( const Expression& cb, const double& j, const double& m, const double& n )
{
  int k_min = std::max(0.,m+n);
  int k_max = std::min(j+m,j+n);
  Expression sum = 0 ; 
  double w2_num = fact(j+m) * fact(j-m) * fact(j+n) * fact(j-n);
  double nc_intpart = 0;
  double ns_intpart = 0; 
  double frac_nc = modf( k_min -(m+n)/2.    , &nc_intpart );
  double frac_ns = modf( j + (m+n)/2. -k_min, &ns_intpart );
  Expression fractional_part = 1 ; 
  if( frac_nc == 0.5 && frac_ns != 0.5 )      fractional_part = fcn::sqrt(1+cb);
  else if( frac_nc != 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb);
  else if( frac_nc == 0.5 && frac_ns == 0.5 ) fractional_part = fcn::sqrt(1-cb*cb);  
  for( double k = k_min; k <= k_max ; ++k )
  {
    double w_den  = fact(k) * fact(j+m-k)*fact(j+n-k) * fact(k-m-n);
    double norm = pow(-1,k) * sqrt(w2_num)/( w_den * pow(2,j) );
    Expression p1 = expandedBinomial( cb, int( k - (m+n)/2.)    );
    Expression p2 = expandedBinomial(-cb, int( j + (m+n)/2. - k)); 
    sum = sum + norm * p1 * p2;  
  }
  auto simplified = NormalOrderedExpression(sum);
  return pow(-1.,j+m) * fractional_part * simplified; 
}

double AmpGen::CG( 
    const double& j1,
    const double& m1,
    const double& j2, 
    const double& m2,
    const double& J,
    const double& M )
{
  if( m1+m2!=M ) return 0;
  double f1 = (2*J+1)*fact(J+j1-j2)*fact(J-j1+j2)*fact(j1+j2-J) ;
  double f2 = fact(j1+m1)*fact(j1-m1)*fact(j2+m2)*fact(j2-m2)*fact(J+M)*fact(J-M);

  double norm = f1 * f2 / fact(J+j1+j2+1) ;
  double sum  = 0;

  for( int nu=0; nu <= j1+j2-J ; ++nu){
    double arg1 = j1+j2-J-double(nu);
    double arg2 = j1  -m1-double(nu);
    double arg3 = j2+  m2-double(nu);
    double arg4 = J-j2+m1+double(nu);
    double arg5 = J-j1-m2+double(nu);
    if( arg1 < 0 || arg2 < 0 || arg3 < 0 || arg4 < 0 || arg5 < 0 ) continue ; 
    int sgn = nu % 2 == 0 ? 1 : -1;
    double to_add =  sgn / (fact(nu)*fact(arg1)*fact(arg2)*fact(arg3)*fact(arg4)*fact(arg5) );
    sum = sum + to_add ; 
  }
  return sqrt(norm) * sum ; 
}

TransformSequence AmpGen::wickTransform( const Tensor& P, 
    const Expression& mass,
    const int& ve,
    const bool& handleZero )
{
  Tensor x({1,0,0}, Tensor::dim(3));
  Tensor y({0,1,0}, Tensor::dim(3));
  Tensor z({0,0,1}, Tensor::dim(3));
  Expression cos_theta = P[2] / fcn::sqrt( P[0]*P[0] + P[1]*P[1] + P[2]*P[2] );
  Expression cos_phi   = P[0] / fcn::sqrt( P[0]*P[0] + P[1]*P[1] );
  Expression sin_phi   = P[1] / fcn::sqrt( P[0]*P[0] + P[1]*P[1] );
    
  Transform rot  = ve == + 1 ? Transform( cos_theta,  sin_phi*x - cos_phi*y, Transform::Type::Rotate) :
                               Transform(-cos_theta, -sin_phi*x + cos_phi*y, Transform::Type::Rotate) ;

  TransformSequence sequence;
  sequence.add(rot);
  if( ve == -1 ) sequence.add( Transform( -1, x, Transform::Type::Rotate ) );
  if( std::real(mass()) != 0. ){
    Transform boost( P[3]/mass, z, Transform::Type::Boost );
    sequence.add(boost); 
  }
  return sequence; 
}

Expression AmpGen::wigner_D(const Tensor& P, 
    const double& J, 
    const double& lA, 
    const double& lB, 
    DebugSymbols* db,
    const std::string& name )
{
  Expression pz = make_cse( P[2] / sqrt( P[0]*P[0] + P[1] * P[1] + P[2]*P[2] ) );  
  Expression pt2 = make_cse( P[0]*P[0] + P[1]*P[1] );
  Expression px = P[0] / sqrt(pt2);
  Expression py = P[1] / sqrt(pt2);

  auto little_d = make_cse ( wigner_d( pz, J, lA, lB ) );
  if( J != 0 && db != nullptr ){
    db->emplace_back("pt2", pt2 );
    db->emplace_back("p2" , sqrt( P[0]*P[0] + P[1] * P[1] + P[2]*P[2] ) );
    db->emplace_back("ϕ("+name+")", atan2( py, px ) );
    db->emplace_back("θ("+name+")", pz );
    db->emplace_back("d[" + std::to_string(J)  +", " + 
                            std::to_string(lA) +", " + 
                            std::to_string(lB) +"](θ)", little_d );
    db->emplace_back("D[" + std::to_string(J)  +", " + 
                            std::to_string(lA) +", " + 
                            std::to_string(lB) +"](Ω)", fpow(px+1i*py,lB-lA) * little_d );
  }
  return  fpow( px + 1i* py, lB - lA ) * little_d; 
}

std::vector<LS> AmpGen::calculate_recoupling_constants( 
    const double& J, 
    const double& M,
    const double& L, 
    const double& S,
    const double& j1,
    const double& j2 ){

  std::vector<LS> rt;
  for( double m1 = -j1; m1 <= j1; ++m1 ){
    for( double m2 = -j2; m2 <= j2; ++m2 ){
      LS f; 
      f.m1 = m1;
      f.m2 = m2;
      f.factor = sqrt( (2.*L + 1. )/( 2.*J + 1. ) );
      f.cg1    = CG(L ,0 ,S ,m1-m2,J,m1-m2);
      f.cg2    = CG(j1,m1,j2,-m2  ,S,m1-m2); 
     // f.cg1    = CG(L ,0 ,S ,m1-m2,J,M);
     // f.cg2    = CG(j1,m1,j2,-m2  ,S,M); 
      f.p      = sqrt( (2*L + 1 )/(2*J+1) );
      f.factor *= f.cg1 * f.cg2;
      if( f.factor != 0 ) rt.push_back(f);
    }
  }
  return rt;
}

Tensor AmpGen::basisSpinor(const int& polState, const int& id)
{
  if(id > 0 && polState ==  1 ) return Tensor({1, 0, 0, 0}, Tensor::dim(4));
  if(id > 0 && polState == -1 ) return Tensor({0, 1, 0, 0}, Tensor::dim(4));
  if(id < 0 && polState ==  1 ) return Tensor({0, 0, 1, 0}, Tensor::dim(4));
  if(id < 0 && polState == -1 ) return Tensor({0, 0, 0, 1}, Tensor::dim(4));
  ERROR("Shouldn't reach here..., polState = " << polState << " id = " << id );
  return Tensor();
}

Tensor AmpGen::basisVector(const int& polState)
{
  double N    = 1./sqrt(2);
  if( polState ==  0 ) return    Tensor(std::vector<complex_t>({0., 0.,1.,0.}), Tensor::dim(4));
  if( polState ==  1 ) return -N*Tensor(std::vector<complex_t>({1., 1i,0.,0.}), Tensor::dim(4));
  if( polState == -1 ) return  N*Tensor(std::vector<complex_t>({1.,-1i,0.,0.}), Tensor::dim(4));
  ERROR("Shouldn't reach here..., polState = " << polState);
  return Tensor();
}

std::vector<LS> userHelicityCouplings( const std::string& key ){
  std::vector<LS> couplings;  
  auto things = NamedParameter<double>( key, 0).getVector();
  if( things.size() % 3 != 0 ) ERROR("Wrong number of tokens");
  for( size_t i = 0 ; i < things.size(); i+=3 ){
    LS coupling; 
    coupling.factor = things[i+0];
    coupling.m1     = things[i+1];
    coupling.m2     = things[i+2];
    couplings.push_back(coupling);
  }
  return couplings;
}


Expression AmpGen::helicityAmplitude(const Particle& particle, 
                                     const TransformSequence& parentFrame, 
                                     const double& Mz, 
                                     DebugSymbols* db, 
                                     int sgn )
{  
  if( particle.daughters().size() > 2 ) return 1; 
  if( particle.daughters().size() == 1 ) 
    return helicityAmplitude( *particle.daughter(0), parentFrame, Mz, db, sgn );
  Tensor::Index a,b,c; 
  auto myFrame = parentFrame; 
  if( particle.spin() == 0 ) myFrame.clear();
  Tensor pInParentFrame = parentFrame(particle.P());
  pInParentFrame.st();
  auto my_sequence = wickTransform(pInParentFrame, fcn::sqrt(particle.massSq()), sgn, true );
  if( ! particle.isHead() ) myFrame.add( my_sequence );
  if( particle.isStable() )
  {
    if( particle.spin() == 0 ) return Mz==0;
    // polarisation spinor / vector etc. in the quantisation of the lab (i.e. along the z-axis or lab particle momentum)
    auto labPol = particle.externalSpinTensor(particle.polState(), db); 
    ADD_DEBUG_TENSOR(labPol, db);
    auto inverseMyTransform = myFrame.inverse();
    if( particle.spin() == 0.5 )
    {
      auto basisSpinor_m1         = basisSpinor( 2*Mz, particle.props()->pdgID() );
      auto labSpinor_m1           = inverseMyTransform( basisSpinor_m1, Transform::Representation::Bispinor );
      ADD_DEBUG_TENSOR(labSpinor_m1, db); 
      return make_cse( Bar(labSpinor_m1)(a)*labPol(a) );
    }
    if( particle.spin() == 1 )
    {
      auto frameVector = basisVector(Mz);
      auto labVector   = inverseMyTransform( frameVector, Transform::Representation::Vector );   
      auto rp = dot( labVector.conjugate(), labPol );
      ADD_DEBUG_TENSOR(labVector, db); 
      ADD_DEBUG(rp, db );
      return rp; 
    }
  }
  auto particle_couplings = particle.spinOrbitCouplings(false);
  auto L = particle.orbital();
  auto& d1 = *particle.daughter(0);
  auto& d2 = *particle.daughter(1);
  double S = 999;
  if( particle.S() == 0 ){ 
    for( auto& l : particle_couplings ) if( l.first == L ){ S = l.second ; break; }
    if( S == 999 ) ERROR("Spin orbital coupling impossible!");
  }
  else S = particle.S()/2.;
  auto recoupling_constants = calculate_recoupling_constants( particle.spin(), Mz, L, S, d1.spin(), d2.spin() );
  auto mod = particle.attribute("helAmp");
  if( mod != stdx::nullopt ) recoupling_constants = userHelicityCouplings( mod.value() );

  if( recoupling_constants.size() == 0 ){    
    WARNING( particle.uniqueString() << " " << particle.spin() << " " << 
        particle.orbitalRange(false).first << " " << particle.orbitalRange(false).second 
        <<  " transition Mz="<< Mz << " to " << d1.spin() << " x " << d2.spin() << " cannot be coupled in (LS) = " << L << ", " << S ); 
    WARNING( "Possible (LS) combinations = " << 
      vectorToString( particle_couplings, ", ", []( auto& ls ){ return "("+std::to_string(int(ls.first)) + ", " + std::to_string(ls.second) +")";} ) );
  } 
  Expression total = 0; 
  for( auto& coupling : recoupling_constants )
  {          
    auto dm = coupling.m1 - coupling.m2;
    if( (d1.name() == "gamma0" && coupling.m1 == 0) || 
        (d2.name() == "gamma0" && coupling.m2 == 0) ) continue;
    auto term = wigner_D(myFrame(d1.P()), particle.spin(), Mz, dm,db, d1.name());
    auto h1   = helicityAmplitude(d1, myFrame, coupling.m1, db, +1);
    auto h2   = helicityAmplitude(d2, myFrame, coupling.m2, db, -1);
    if( db != nullptr ){ 
      db->emplace_back( "coupling" , coupling.factor );
      if( coupling.factor != 1 ) db->emplace_back( "C x DD'", coupling.factor * term * h1 * h2 );
    }
    total = total + coupling.factor * term * h1 * h2;
  }
  return total;
}
