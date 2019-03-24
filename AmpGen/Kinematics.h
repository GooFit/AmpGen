#ifndef AMPGEN_KINEMATICS_H
#define AMPGEN_KINEMATICS_H
#include <stddef.h>
#include <array>
#include <tuple>
#include <vector>

#include "TLorentzVector.h"
#include "TVector3.h"

namespace AmpGen
{
  class Event; 

  /** @defgroup Kin Kinematics
      @brief Functionality related to kinematics.
      Assorted functors for calculating kinematic quantities on events, such as helicity cosines and acoplanarities. 
      Also contains utilities for boosting and rotating four vectors and building general transform sequences on different Lorentz objects. 
   */

   /** @ingroup Kin class HelicityCosine 
      @brief Functor to compute the angle between set of particles {1} and {2} in the rest frame of set {3}. 
    */
  class HelicityCosine {
    public:
      HelicityCosine(const std::vector<size_t>& p1, const std::vector<size_t>& p2,
          const std::vector<size_t>& pR);

      HelicityCosine(const size_t& i, const size_t& j, const std::vector<size_t>& pR);

      double operator()( std::vector<Event>::iterator evt ) const;
      double operator()( const Event& evt ) const;
    private:
      std::vector<size_t> _i, _j, _pR;
  };

  /** @ingroup Kin class MomentumTransfer
      @brief Functor to calculate the linear momemtum between particles {1} and {2} in the rest frame of {1} + {2}. 
   */ 
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

  /** @ingroup Kin function acoplanarity
      @brief The extent to which a four body decay occurs within a single decay frame. 
      Defined by the angle between the normals of decay planes of two of the quasi two-body subsystems, i.e. 
      @f[
        \chi = \cos^{-1} \left( \frac{ ( p_1 \times p_2 ) \cdot ( p_3 \times p_4 ) }{ | p_1 \times p_2 | | p_3 \times p_4 | } \right),
      @f]
     where each of the three-vectors is calculated in the rest frame of the decaying particle. 
    */
  double acoplanarity( const Event& evt );

  double PHI( const Event& evt );
  double phi( const Event& evt, int i, int j, int k, int w );

  std::vector<double> rotate( const std::vector<double>& input, const std::vector<double>& n, const double& v );
  void boost( Event& evt, const std::tuple<double, double, double>& n, const double& v );
  void rotate( Event& evt, const std::tuple<double, double, double>& n, const double& v );

  void rotateBasis( Event& evt, const TVector3& p1, const TVector3& p2, const TVector3& p3 );
  
  /** @ingroup Kin function dotProduct 
      @brief Helper function to calculate the (space-like) dot product between vectors p1 and p2 in the rest frame of pX.  
      Helper function to calculate the (space-like) dot product between vectors p1 and p2 in the rest frame of pX, which is given by
      @f[
         d = \left( -g_{\mu\nu} + \frac{p_X^{\mu}p_X^{\nu}}{p_X^2} \right) p_1^{\mu} p_2^{\nu}
      @f]
    */
  double dotProduct( const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& pX );
  
  /** @ingroup Kin function pFromEvent
      @brief Helper function to extract a TLorentzVector of the ith decay product from an AmpGen::Event.
    */
  TLorentzVector pFromEvent( const Event& evt, const size_t& ref );

  /** @ingroup Kin function pFromEvent
      @brief Helper function to extract a TLorentzVector of the sum of {2} decay product from an AmpGen::Event.
    */
  TLorentzVector pFromEvent( const Event& evt, const std::vector<size_t>& ref );
  
} // namespace AmpGen
#endif
