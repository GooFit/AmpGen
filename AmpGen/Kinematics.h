#ifndef AMPGEN_KINEMATICS_H
#define AMPGEN_KINEMATICS_H
#include <array>
#include <tuple>
#include <vector>

#include "AmpGen/EventList.h"
#include "TLorentzVector.h"
namespace AmpGen
{
  /// \defgroup Kin Kinematics
  /// Assorted functors for computing kinematic quantities on events, such as helicity cosines and acoplanarities. 
  /// Also contains utilities for boosting and rotating events  

  /// \ingroup Kin class HelicityCosine 
  /// \brief Functor to compute the angle between set of particles {1} and {2} in the rest frame of set {3}. 
//  enum RepresentationType {
//    scalar, vector, spinor, antispinor, tensor; 
//  };

  class HelicityCosine {
    public:
      HelicityCosine( const std::vector<size_t>& p1, const std::vector<size_t>& p2,
          const std::vector<size_t>& pR );

      HelicityCosine( const size_t& i, const size_t& j, const std::vector<size_t>& pR );

      double operator()( std::vector<Event>::iterator evt ) const;
      double operator()( const Event& evt ) const;
    private:
      std::vector<size_t> _i, _j, _pR;
  };

  class MomentumTransfer {
    public:
      MomentumTransfer( const std::vector<size_t>& _p1, const std::vector<size_t>& _p2 );
      double operator()( const Event& evt ) const;
    private: 
      double Q2( const double& s, const double& s1, const double& s2 ) const;
      std::vector<size_t> p1;
      std::vector<size_t> p2;
      std::vector<size_t> s;
  };

  TLorentzVector pFromEvent( const Event& evt, const size_t& ref );
  TLorentzVector pFromEvent( const Event& evt, const std::vector<size_t>& ref );

  double acoplanarity( const Event& evt );
  double Product( const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& pX );

  double trihedralAngle( const Event& evt );
  double TripleProduct( const Event& evt );

  double PHI( const Event& evt );
  double phi( const Event& evt, int i, int j, int k, int w );

  std::vector<double> rotate( const std::vector<double>& input, const std::vector<double>& n, const double& v );
  void boost( Event& evt, const std::tuple<double, double, double>& n, const double& v );
  void rotate( Event& evt, const std::tuple<double, double, double>& n, const double& v );

  void rotateBasis( Event& evt, const TVector3& p1, const TVector3& p2, const TVector3& p3 );
  
  Tensor BoostMatrix( const Tensor& p );
} // namespace AmpGen
#endif
