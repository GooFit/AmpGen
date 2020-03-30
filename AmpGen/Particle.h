#ifndef AMPGEN_PARTICLE_H
#define AMPGEN_PARTICLE_H

#include <complex>
#include <vector>
#include <memory>

// hack to include optional from https://codereview.stackexchange.com/questions/136350/seamlessly-migrating-experimental-optional-to-optional
#if __cplusplus >= 201703L
#include <optional>
namespace stdx {
  using namespace ::std;
}
#elif __cplusplus >= 201402L
  #include  <experimental/optional>
  namespace stdx {
    using namespace ::std;
    using namespace ::std::experimental;
  }
#else 
#   error "Require c++ std >=14"
#endif


#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/QuarkContent.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/enum.h"
namespace AmpGen
{
  /** @class Particle
      @brief Describes a particle, its decay process and subsequent decay products, which are also Particles. 

    Encodes a multi-body decay tree structure, is largely limited to describe a sequence 
    of quasi two-body processes (i.e. the isobar model) by the implemented lineshapes (propagators) and 
    spin structure (vertices) 
    Decay chains are usually constucted by parsing strings, for example:

    \code{cpp}
    D0{rho(770)0{pi+,pi-},K0S0}
    \endcode

    is the tree structure for the decay of a neutral D-meson into a @f$\rho@f$ meson and the short-lived neutral kaon, 
    @f$ K^0_S @f$. States are lists of particles, separated by commas, and included within {} braces. 
    As a trivial example, this decay descriptor can be used to generate some of its own documentation, in this case a formatted latex string:

    \code{cpp}
    AmpGen::Particle example("D0{rho(770)0{pi+,pi-},K0S0}");
    std::cout << example.texLabel(true) << std::endl;
    \endcode

    produces the latex source for the output 

    @f[
    D^{0}\rightarrow \rho(770)^{0}\left[\pi^{+} \pi^{-} \right] K_{S}^{0}
    @f]

    There are also modifiers to the amplitude that alter the details of the decay process. 
    For example, there may multiple orbital angular momentum substates. 
    These will be denoted by the use of [] braces, so for example 

    \code{cpp}
    a(1)(1260)+{rho(770)0,pi+}
    a(1)(1260)+[D]{rho(770)0,pi+} 
    \endcode

    correspond to the S-wave and D-wave decays of the @f$a(1)(1260)^{+}@f$ meson.
    As a rule, the lowest orbital state permitted by the relevant conservation laws of the decay is 
    used if the orbital state is not specified, so the conservation of angular momentum, 
    and the conservation of parity if the decay proceeds via the strong or electromagnetic force.

    The modifier syntax is also used to specify a different choice of lineshape for the resonance. 
    For example, a common parameterisation for the @f$\rho(770)@f$ meson is the Gounaris-Sakurai propagator, 
    which accounts for dispersive corrections to the @f$I=1@f$ dipion system. 
    In this example

    \code{cpp}
    rho(770)0[GounarisSakurai]{pi+,pi-}
    \endcode 

    Multiple modifiers can be applied to the same particle, by including them in a semi-colon separated list. 
    For example, for three-body decays of broad resonances, such as the @f$a_1(1260)@f$, the propagator must take the evolution of 
    the width of resonance from a numerical calculation, which can be supplied via a cubic spline (hence using the GSpline propagator).
    Additionally, for the decay chain @f$ a_1(1260)^{+} \to \rho^{0} \pi^{+} @f$, the decay products can either be in a relative S or D wave, and hence 
    the two particle descriptors

    \code{cpp}
    a(1)(1260)+[GSpline]{rho(770)0,pi+}
    a(1)(1260)+[D;GSpline]{rho(770)0,pi+}
    \endcode

    are relevant for the decay.   
    Similar to other components of AmpGen, Particles will rarely be constructed in the C++ context, 
    and will instead be instantiated dynamically at runtime from a user supplied options file. */
  class ParticleProperties; 
      
  declare_enum( spinFormalism, Covariant, Canonical )
  declare_enum( spinBasis    , Dirac    , Weyl ) 

  class Particle
  {
    public:
      /// Default Constructor
      Particle();
      
      /// Constructor that takes a pair of other particles (i.e. this particle's decay products) as arguments and looks up the properties of this particle using the particle name. 
      Particle( const std::string& name, const Particle& p1, const Particle& p2 ); 
      
      /// Constructor that takes a pair of other particles (i.e. this particle's decay products) as arguments and looks up the properties of this particle using the PDG MC ID.       
      Particle( const int& pdg_id, const Particle& p1, const Particle& p2 ); 
      
      /// Constructor by name and with an index to match to the event type
      Particle( const std::string& name, const unsigned int& index ); 

      /// Constructor that takes a decay descriptor as an argument and a list of final state particles to match to the event type. Constructs the entire decay tree.  
      Particle( const std::string& decayString, const std::vector<std::string>& finalStates = {}, const bool& orderDaughters = true );

      /// (Quasi) Constructor that returns the (quasi)CP conjugated amplitude. The full behaviour of the amplitude is made more complicated by the ordering convention. 
      Particle conj(bool invertHead = true, bool reorder = true);

      /// CP conjugate this particle //

      void conjThis();

      static bool isValidDecayDescriptor( const std::string& decayDescriptor ); 
      /// Set the orbital quantum number 'L' for this decay.       
      void setOrbital( const unsigned int& orbital );

      /// Set the lineshape for the decay of this particle. 
      void setLineshape( const std::string& lineshape );

      /// Set the index'th daughter of this to particle. 
      void setDaughter( const Particle& particle, const unsigned int& index );
 
      /// Set the parent particle for this particle 
      void setParent( const Particle* particle );

      /// Set the index of this particle, i.e. where it is positioned in the event data structure. 
      void setIndex( const unsigned int& index, const bool& setOri = false );

      /// Remove all of the decay products of this particle
      void clearDecayProducts();

      /// Add some modifier to the particle, such as a lineshape or a different spin state
      void addModifier( const std::string& mod );

      /// Parse some set of modifiers, delimited with semicolons.
      void parseModifier( const std::string& mod );

      /// Set some particle ordering of the decay products of this particle, mostly used internally by the symmetrisation
      void setOrdering( const std::vector<size_t>& ordering );

      /// Set the particle name 
      void setName(const std::string& name);

      /// Add a decay product
      void addDaughter( const std::shared_ptr<Particle>& particle );

      /// Set the polarisation state of this particle, which is twice the projection of the spin along the quantisation axis.
      void setPolarisationState( const int& state );
      void setPolarisationState( const std::vector<int>& state);

      /// Returns the range of orbital angular momentum between the decay products
      std::pair<size_t,size_t> orbitalRange( const bool& converseParity = true ) const;

      /// Returns the set of possible spin-orbit couplings allowed by conservation of angular momentum, and if specified parity
      std::vector<std::pair<double,double>> spinOrbitCouplings( const bool& conserveParity = true ) const;

      /// Return the additional optional attribute keyed by variable key 
      stdx::optional<std::string> attribute(const std::string& key) const; 

      /// Return the particleProperties object for this particle
      const ParticleProperties* props() const;
      
      QuarkContent quarks() const;          ///< Return the quarks of this particle

      QuarkContent daughterQuarks() const; ///< Returns the quark content of the sum of the decay products of this particle

      int parity() const;            ///< Returns the parity of this particle
      int finalStateParity() const;  ///< Returns the parity of the final state of this particle 

      int polState() const;          ///< Returns the polarisation state, i.e. twice the projection of the spin along the quantisation axis, of this particle. 
      int CP() const;                ///< Returns the CP of this decay. 
      int C() const;                 ///< Returns the C quantum number for this decay
      double mass() const;           ///< Returns the (PDG) mass of the particle
      double spin() const;           ///< Returns the spin of the particle
      double S() const;              ///< Returns the spin configuration of the decay products of the particle
      unsigned L() const;            ///< Returns the orbital angular 

      bool isHead() const;           ///< Returns whether if this particle is the head of the decay, i.e. has no parent
      bool isWeakDecay() const;      ///< Returns whether is this particle decays weakly or not
      bool isStateGood() const;      ///< Returns whether this particle, and its decays have been configured correctly 
      bool isStable() const;         ///< Check whether this particle is stable, has any decay products. 
      bool isQuasiStable() const;    ///< Check whether the particle is quasi-stable, i.e. may have some appreciable flight distance
      bool conservesParity( unsigned int L = 0 ) const; ///< Check whether the decay of this particle with angular momentum L conserves parity or not

      unsigned index() const;             ///< Returns the current index of the particle in event data structure. Can differ from the original index due to symmetrisation
      unsigned originalIndex() const;     ///< Returns the original index of the particle
      std::string name() const;           ///< Name of the decaying particle.
      std::string lineshape() const;      ///< Name of the propagator to use for the decay of this particle.

      std::string uniqueString() const;   ///< Returns the unique string (i.e. decay descriptor) that identifies this decay, which can be parsed to generate the decay tree.
      std::string decayDescriptor() const;///< Returns the unique string (i.e. decay descriptor) that identifies this decay, which can be parsed to generate the decay tree.
      
      /// The string that describes the spin/orbital topology of this decay, 
      /// i.e. replacing specific particle names with their spins.
      std::string topologicalString() const;
      
      /// The string that describes the spin/orbit configuration of this decay.
      std::string orbitalString() const;
      
      /// Decay descriptor formatted as LaTeX for this decay. 
      std::string texLabel( const bool& printHead = false, const bool& recurse=true ) const;
      
      /// Return the eventType for this decay (i.e. the initial and final state particles) 
      EventType eventType() const;

      /// Returns the indexth decay product of this particle
      std::shared_ptr<Particle> daughter( const size_t& index );
      
      /// Returns in indexth decay product of this particle (as constant)
      std::shared_ptr<Particle> daughter( const size_t& index ) const;
      
      /// Returns in indexth decay product of this particle (as constant)
      std::shared_ptr<Particle> daughter(const std::string& name, const int& maxDepth=-1) const;
     
      /// Vector of decay products of this particle 
      std::vector<std::shared_ptr<Particle>> daughters() const;

      /// Get orderings of the final state that are identical to each other, i.e. those that only differ by exchanging identical particles. 
      std::vector<std::vector<size_t>> identicalDaughterOrderings() const;
      
      /// Returns the final state particles for this decay process.
      std::vector<std::shared_ptr<Particle>> getFinalStateParticles( const bool& sort = true ) const;
      
      /// Calculate the particle tree only including quasi-stable processes, 
      /// i.e. only includes states with lifetimes > ParticleProperties::qsThreshold (default ~ 1 KeV ) 
      Particle quasiStableTree() const;

      /// Calculates the momentum sum of the decay products
      Tensor P() const;
    
      /// Calculates the momentum difference between the decay products (only well defined for quasi two-body processes )
      Tensor Q() const;
     
      /// Calculates the spin tensor or generalised current for this particle 
      Tensor spinTensor( DebugSymbols* db = nullptr ) const;
      
      /// Calculates the polarisation vector / spinor etc. of this particle, used for the initial/final state particles 
      Tensor externalSpinTensor(const int& polState, DebugSymbols* db = nullptr) const;
      
      /// Calculates the invariant mass-squared of the mass of this particle
      Expression massSq() const;                                

      /// Calculates the lineshape / propagator for this particle. 
      Expression propagator( DebugSymbols* db = nullptr ) const;

      /// Calculates the total expression for this particle, including symmetrisation and the current polarisation state
      Expression getExpression( DebugSymbols* db = nullptr, const unsigned int& index = 0 );


      /// Experimental Zemach Code
      
      //Factorial function - helper for Zemach
      int Factorial(int n){
      	if (n==0){
      		return 1;
      	}
      	else{
      		return n * Factorial(n-1);
      	}
      }
      
      Expression LegendreZemach(Expression x, double n){
      
      	Expression L = 0;
      	if (n==0){
      		L = 1;
      	}
      	else if (n==1){
      		L = x;
      	}
      	else{
      		L = (2 * n -1) * LegendreZemach(x, n-1) - (n-1)*LegendreZemach(x, n-2);
      		L = L/n;
      	}
      	return  L;
      
      }
      
      Expression cosHel(const std::shared_ptr<Particle>& R, const std::shared_ptr<Particle>& h){
      
      	auto pR = R->P();
      	auto ph = h->P();
      	auto R1 = R->daughters()[0];
      	auto R2 = R->daughters()[1];
      	auto pR1 = R1->P();
      	auto pR2 = R2->P();
      
      	auto pD = pR + ph;
      
      	auto m2D = dot(pD, pD);
      	auto m2R = dot(pR, pR);
      	auto m2R1 = dot(pR1, pR1);
      	auto m2R2 = dot(pR2, pR2);
      	auto m2h = dot(ph, ph);
      	auto m2R1h = dot(pR1 + ph, pR1 + ph);
      
      	auto EhR = (m2D - m2R - m2h)/(2 * fcn::sqrt(m2R));
      	auto ER1R = (m2R + m2R1 - m2R2)/(2 * fcn::sqrt(m2R));
      
      	auto PR1R = fcn::sqrt(ER1R * ER1R - m2R1);
      	auto PhR = fcn::sqrt(EhR * EhR - m2h);
      
      	auto cosH = (m2R1h - m2R1 - m2h - 2 * ER1R * EhR)/(2 * PR1R * PhR);
      	return cosH;
      
      }
      
      Expression pq( const std::shared_ptr<Particle>& R, const std::shared_ptr<Particle>& h, bool parent){
      	auto pR = R->P();
      	auto ph = h->P();
      	auto R1 = R->daughters()[0];
      	auto R2 = R->daughters()[1];
      	auto pR1 = R1->P();
      	auto pR2 = R2->P();
      
      	auto pD = pR + ph;
      
      	auto m2D = dot(pD, pD);
      	auto m2R = dot(pR, pR);
      	auto m2R1 = dot(pR1, pR1);
      	auto m2R2 = dot(pR2, pR2);
      	auto m2h = dot(ph, ph);
      	auto m2R1h = dot(pR1 + ph, pR1 + ph);
      
      	auto EhR = (m2D - m2R - m2h)/(2 * fcn::sqrt(m2R));
      	auto ER1R = (m2R + m2R1 - m2R2)/(2 * fcn::sqrt(m2R));
      
      	auto PR1R = fcn::sqrt(ER1R * ER1R - m2R1);
      	auto PhR = fcn::sqrt(EhR * EhR - m2h);
      
      	auto p = PhR;
      	auto q = PR1R;
      
      	auto EhD = (m2D + m2h - m2R)/(2 * fcn::sqrt(m2D));
      	auto PhD = fcn::sqrt(EhD * EhD - m2h);
      	auto pstar = PhD;
      	if (parent){
      		return pstar*q;
      	}
      	else {
      		return p*q;
      	}
      
      
      }
      
          //Zemach Tensor stolen from A. Poluektov's TFA package - see ZemachTensor in TensorFlowAnalysis/Kinematics.py
          //D -> ABC at the moment
      Expression Zemach(Expression m2ab, Expression m2ac, Expression m2bc, Expression m2d, Expression m2a, Expression m2b, Expression m2c, double spin){
      	Expression Z = 1;
      	if (spin==0){
      		Z = 1;
      	}
      	else if(spin==1){
      		Z = m2ac-m2bc+(m2d-m2c)*(m2b-m2a)/m2ab;
      	}
      	else if (spin ==2){
      		Z = (m2bc-m2ac+(m2d-m2c)*((m2a-m2b)/m2ab) * ((m2a-m2b)/m2ab) -1./3.*(m2ab-2.*(m2d+m2c) + (m2d-m2c) * (m2d-m2c)/m2ab)*(m2ab-2.*(m2a+m2b)+(m2a-m2b)*(m2a - m2b)/m2ab));
      	}
      	return Z;
      	
      }
      
      Expression ZemachLaura(const std::shared_ptr<Particle>& R, const std::shared_ptr<Particle>& h, double spin, bool parent){
      	auto pq_ = pq(R, h, parent);
      	auto cosH = cosHel(R, h);
      
      	Expression Z=1;
      	if (spin == 0){
      		Z = 1;
      	}
      	else {
      		auto c = pow(-2, spin) * Factorial(spin)/Factorial(Factorial(2 * spin - 1));
      		Z *= c;
      		for (int i=0; i < spin; i++){
      			Z *= pq_;
      		}
      		Z*=LegendreZemach(cosH, spin);
      
      
      	}
      	return Z;
      }
      
      Expression ZemachGooFit(const std::shared_ptr<Particle>& R, const std::shared_ptr<Particle>& C, double spin){
      	if (spin==0){
      		return 1;
      	}
      	auto A = R->daughters()[0];
      	auto B = R->daughters()[1];
      
      	auto pA = A->P();
      	auto pB = B->P();
      	auto pC = C->P();
      	auto pD = pA + pB + pC;
      
      	Expression mAA = dot(pA + pA, pA + pA);
      	Expression mBB = dot(pB + pB, pB + pB);
      	Expression mCC = dot(pC + pC, pC + pC);
      	Expression mDD = dot(pD + pD, pD + pD);
      
      	Expression mAB = dot(pA + pB, pA + pB);
      	Expression mBC = dot(pB + pC, pB + pC);
      	Expression mAC = dot(pA + pC, pA + pC);
      
      	//Mass Factor
      	Expression MF = 1/mAB;
      	
      	//Spin=1 - also need for spin 2
      	
      	Expression spinFactor1 = -1;
      	spinFactor1 = ((mBC - mAC) + MF * (mDD - mCC) * (mAA - mBB)) * spinFactor1;
      
      	if (spin==1){
      
      		return spinFactor1;
      	}
      
      	auto spinFactor2 = spinFactor1 * spinFactor1;
      
      	auto extraTerm = mAB - 2 * mDD - 2*mCC + MF * (mDD - mCC) * (mDD - mCC);
      	extraTerm = (mAB - 2*mAA - 2*mBB + MF * (mAA - mBB) * (mAA - mBB)) * extraTerm/3.;
      
      	spinFactor2 = spinFactor2 - extraTerm;
      
      	if (spin==2){
      		return spinFactor2;
      	}
      
      	return 1;
      
      
      	
      
      }
      
      
      Expression GooFitCMom(Expression m2res, Expression m1, Expression m2){
      	Expression k1 = 1 - fcn::pow(m1 + m2, 2)/m2res;
      	Expression k2 = 1 - fcn::pow(m1 - m2, 2)/m2res;
      	Expression cMom = 1/2 * fcn::sqrt(m2res * k1 * k2);
      	return cMom;
      }
      
      Expression GooFitDamping(const std::shared_ptr<Particle>& R, const std::shared_ptr<Particle>& C, double spin){
      	auto A = R->daughters()[0];
      	auto B = R->daughters()[1];
      
      	auto pA = A->P();
      	auto pB = B->P();
      	auto pC = C->P();
      	auto pD = pA + pB + pC;
      
      	Expression mAA = dot(pA + pA, pA + pA);
      	Expression mBB = dot(pB + pB, pB + pB);
      	Expression mCC = dot(pC + pC, pC + pC);
      	Expression mDD = dot(pD + pD, pD + pD);
      
      	Expression mAB = dot(pA + pB, pA + pB);
      	Expression mBC = dot(pB + pC, pB + pC);
      	Expression mAC = dot(pA + pC, pA + pC);
      
      	Expression mA = fcn::sqrt(mAA);
      	Expression mB = fcn::sqrt(mBB);
      	Expression mC = fcn::sqrt(mCC);
      	Expression mD = fcn::sqrt(mDD);
      
      	Expression rR = 5;
      
      	Expression Cmom1 = GooFitCMom(R->massSq(), mA, mB);
      	Expression Cmom2 = GooFitCMom(mAB, mA, mB);
      
      	Expression sq1 = fcn::pow(rR*Cmom1, 2);
      	Expression sq2 = fcn::pow(rR*Cmom2, 2);
      
      	Expression df1 = 1 + sq1;
      	Expression df2 = 1 + sq2;
      
      	if (spin == 2){
      		df1 = df1 + 8 + 2 * sq1 * sq1;
      		df2 = df2 + 8 + 2 * sq2 * sq2;
      	}
      	return df2/df1;
      
      
      	
      
      }




      /// Calculate the transition matrix for this decay 
      Tensor transitionMatrix( DebugSymbols* db = nullptr );
      bool operator<( const Particle& other );
      bool operator>( const Particle& other );

      enum MatchState 
      {
        None                  = ( 1<<0 ), 
        Exact                 = ( 1<<1 ), 
        PartialExpansion      = ( 1<<2 ), 
        DifferentOrbital      = ( 1<<3 ), 
        DifferentPolarisation = ( 1<<4 )
      };
      /// matches Check the matching between two decay chains, according to the MatchState enum. 

      unsigned int matches( const Particle& other ) const; 
      std::string makeUniqueString();                        ///< Generate the decay descriptor for this decay. 
    private:
      std::string m_name                     = {""};         ///< Name of the particle
      const ParticleProperties* m_props      = {nullptr};    ///< Particle Properties from the PDG
      std::string m_lineshape                = {"BW"};       ///< Propagator to use
      std::string m_uniqueString             = {""};         ///< Unique string of particle tree
      int m_parity                           = {0};          ///< Intrinsic parity of particle
      int m_polState                         = {0};          ///< Projection of the spin along the quantisation axis, i.e. 'z'
      unsigned m_index                       = {999};        ///< Index, for constructing four-momenta
      unsigned m_originalIndex               = {999};        ///< Starting index, used in Bose-symmetrisation
      unsigned m_orbital                     = {0};          ///< Orbital angular momentum between daughters
      unsigned m_spinConfigurationNumber     = {0};          ///< Spin configuration quantum number 'S'
      unsigned m_minL                        = {0};          ///< Minimum orbital angular momentum
      bool m_usesDefaultLineshape            = {false};      ///< Flag to check if default shape is used
      bool m_isStateGood                     = {true};       ///< Flag to check the decay is well-formed
      std::vector<std::shared_ptr<Particle>> m_daughters;    ///< Array of daughter particles
      std::vector<std::string> m_modifiers;                  ///< Additional modifiers for amplitude
      const Particle*   m_parent             = {nullptr};    ///< Pointer to the parent particle of this particle
      void pdgLookup();                                      ///< Lookup information from the PDG database (using ParticlePropertiesList)
      bool hasModifier( const std::string& modifier ) const; ///< Check if this particle has a given modifier
      std::string modifierString() const;                    ///< Re-generate modifier string used to create particle
      void sortDaughters();                                  ///< Recursively order the particle's decay products. 

      NamedParameter<spinFormalism> m_spinFormalism  = {"Particle::SpinFormalism"  ,spinFormalism::Covariant, optionalHelpString("Formalism to use for spin calculations", {  
             {"Covariant", "[default] Covariant Tensor, based on Rarita-Schwinger constraints on the allowed covariant wavefunctions."}
           , {"Canonical", "Canonical formulation, based on rotational properties of wavefunctions, i.e. Wigner D-matrices and Clebsch-Gordan for (L,S) expansion."} } ) };

      NamedParameter<spinBasis>     m_spinBasis      = {"Particle::SpinBasis", spinBasis::Dirac, optionalHelpString("Basis to use for calculating external polarisation tensors / spinors.", {
                      {"Dirac", "[default] Quantises along the z-axis"}
                    , {"Weyl" , "Quantises along the direction of motion"}} )};
      NamedParameter<std::string> m_defaultModifier = {"Particle::DefaultModifier","", "Default modifier to use for lineshapes, for example to use normalised vs unnormalised Blatt-Weisskopf factors."};
  };
  std::ostream& operator<<( std::ostream& os, const Particle& particle );
} // namespace AmpGen

#endif
