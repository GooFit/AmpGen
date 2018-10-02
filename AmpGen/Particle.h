#ifndef AMPGEN_PARTICLE_H
#define AMPGEN_PARTICLE_H

/// STL
#include <complex>
#include <vector>

/// AmpGen
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Tensor.h"

namespace AmpGen
{
  /** @class Particle
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
        Delta(1232)++{p+,pi+}
        Delta(1232)++[D]{p+,pi+}
      \endcode

      correspond to the S-wave and D-wave decays of the @f$\Delta(1232)^{++}@f$ baryon.
      As a rule, the lowest orbital state permitted by the relevant conservation laws of the decay is 
      used if the orbital state is not specified, so the conservation of angular momentum, 
      and the conservation of parity if the decay proceeds via the strong or electromagnetic force.
   
      The modifier syntax is also used to specify a different choice of lineshape for the resonance. 
      For example, a common parameterisation for the @f$\rho(770)@f$ meson is the Gounaris-Sakurai propagator, 
      which accounts for dispersive corrections to the @f$I=1@f$ dipion scattering. In this example

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
      and will instead be instantiated dynamically at runtime from a user supplied options file. 
   */
  class Particle
  {
  private:
    const ParticleProperties* m_props;                     ///< Particle Properties from the PDG
    std::string m_name;                                    ///< Name of the particle
    std::string m_lineshape;                               ///< Propagator to use
    std::string m_uniqueString;                            ///< Unique string of particle tree
    int m_parity;                                          ///< Intrinsic parity of particle
    int m_polState;                                        ///< polarisation state 
    unsigned int m_index;                                  ///< Index, for constructing four-momenta
    unsigned int m_originalIndex;                          ///< Starting index, used in Bose-symmetrisation
    unsigned int m_orbital;                                ///< Orbital angular momentum between daughters
    unsigned int m_spinConfigurationNumber;                ///< Spin configuration quantum number
    unsigned int m_minL;                                   ///< Minimum orbital angular momentum
    bool m_isHead;                                         ///< Flag that particle is head of decay chain
    bool m_isStateGood;                                    ///< Flag to check the decay is well-formed
    bool m_usesDefaultLineshape;                           ///< Flag to check if default shape is used
    std::vector<std::shared_ptr<Particle>> m_daughters;    ///< Array of daughter particles
    std::vector<std::string> m_modifiers;                  ///< Additional modifiers for amplitude

    void pdgLookup();                                      ///< Lookup information from the PDG database (using ParticlePropertiesList)
    bool hasModifier( const std::string& modifier ) const; ///< Check if this particle has a given modifier
    std::string modifierString() const;                    ///< Re-generate modifier string used to create particle
    void sortDaughters();                                  ///< Recursively order the particle's decay products. 
  public:
    Particle();                                            
    Particle( const std::string& name, const Particle& p1, const Particle& p2 ); 
      ///< Constructor that takes a pair of other particles (i.e. this particle's decay products) as arguments and looks up the properties of this particle using the particle name. 
    Particle( const int& pdg_id, const Particle& p1, const Particle& p2 ); 
      ///< Constructor that takes a pair of other particles (i.e. this particle's decay products) as arguments and looks up the properties of this particle using the PDG MC ID.       
    Particle( const std::string& name, const unsigned int& index );        
      ///<Constructor by name and with an index to match to the event type
    Particle( const std::string& decayString, const std::vector<std::string>& finalStates = {}, const bool& orderDaughters = true );
      ///<Constructor that takes a decayString as an argument and a list of final state particles to match to the event type. Constructs the entire decay tree.  

    void setOrbital( const unsigned int& orbital );
    void setLineshape( const std::string& lineshape );
    void setDaughter( const Particle& particle, const unsigned int& index );
    void setTop( bool state = true );
    void setIndex( const unsigned int& index, const bool& setOri = false );

    void addModifier( const std::string& mod );
    void parseModifier( const std::string& mod );
    int conjugate( bool invertHead = false , bool reorder = true);
    void setOrdering( const std::vector<size_t>& ordering );
    void addDaughter( const std::shared_ptr<Particle>& particle );
    void setPolarisationState( const int& state );
    std::pair<size_t,size_t> orbitalRange( bool converseParity = true ) const; ///< Range of possible orbital angular momenta between decay products

    const ParticleProperties* props() const;
    MultiQuarkContent quarks() const;
    MultiQuarkContent daughterQuarks() const;
    int parity() const;
    int finalStateParity() const;
    int polState() const;
    double mass() const;
    double spin() const;    
    bool isTop() const;
    bool isWeakDecay() const;
    bool isStateGood() const;
    bool isStable() const;
    bool isQuasiStable() const;
    bool conservesParity( unsigned int L = 0 ) const;
    bool checkExists() const;

    unsigned int orbital() const;
    unsigned int index() const;
    unsigned int originalIndex() const;

    std::string name() const;                                            ///< Name of the decaying particle.
    std::string lineshape() const;                                       ///< Name of the propagator to use for the decay of this particle.
    std::string vertexName() const;
    
    std::string uniqueString() const;                                    ///< Returns the unique string that identifies this decay / can be parsed to generate the decay tree.
    std::string topologicalString() const;                               ///< String that describes the angular momentum configuration. 
    std::string orbitalString() const;                                   ///< String that describes the orbital configuration of this decay (S,P,D) + possible modifiers. 
    std::string texLabel( const bool& printHead = false ) const;         ///< Get the TeX formatted label for this decay process.
    std::string makeUniqueString();                                      ///< Generate the unique string that identifies this decay / can be parsed into this structure.
    EventType eventType() const;                                         ///< EventType 

    std::shared_ptr<Particle> daughter( const size_t& index );           ///< Return the ith decay product
    std::shared_ptr<Particle> daughter( const size_t& index ) const;     ///< Return the ith decay product
    std::vector<std::shared_ptr<Particle>> daughters() const;            ///< vector of decay products

    std::vector<std::vector<size_t>> identicalDaughterOrderings() const; ///< Get orderings of the final state that are identical to each other, i.e. those that only differ by exchanging identical particles. 
    std::vector<std::shared_ptr<Particle>>
    getFinalStateParticles( const bool& sort = true ) const; ///< Retrieve vector of final state particles.
    Particle quasiStableTree() const;                        ///< Returns the tree of quasi-stable processes, i.e. only includes states with lifetimes > ParticleProperties::qsThreshold (default ~ 1 KeV ) 
    Tensor P() const;                                      ///< momentum sum of daughters
    Tensor Q() const;                                      ///< momentum difference between daughters (only well-defined for quasi two-body processes)
    Tensor SpinTensor( DebugSymbols* db = nullptr ) const; ///< Spin current associated with decay (including currents
                                                           ///< and matrix elements of decay products)
    Tensor ExternalSpinTensor(const int& polState, DebugSymbols* db = nullptr) const;   ///< Spin tensor associated with an external particle that carries spin, i.e. fermions and photons.
    Expression massSq() const;                                ///< Invariant mass-squared of particle
    Expression Lineshape( DebugSymbols* db = nullptr ) const; ///< Lineshape expression of particle
    Expression getExpression( DebugSymbols* db = nullptr, const unsigned int& index = 0 );

    bool operator<( const Particle& other );

    enum MatchState 
    {
      None                  = ( 1<<0 ), 
      Exact                 = ( 1<<1 ), 
      PartialExpansion      = ( 1<<2 ), 
      DifferentOrbital      = ( 1<<3 ), 
      DifferentPolarisation = ( 1<<4 )
    };
    unsigned int matches( const Particle& other ) const; ///< Check the matching between two decay chains, according to the MatchState enum. 
  };
} // namespace AmpGen

#endif
