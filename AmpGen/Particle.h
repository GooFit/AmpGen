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
  ///\class Particle
  ///Encodes a multi-body decay tree structure, essentially only limited by
  ///spin factor and propagators implemented to quasi two-body processes,
  ///and is the central object for computing the expression tree of
  ///a matrix element
  
  class Particle
  {
  private:
    const ParticleProperties* m_props;                  ///< Particle Properties from the PDG
    std::string m_name;                                 ///< Name of the particle
    std::string m_lineshape;                            ///< Propagator to use
    std::string m_uniqueString;                         ///< Unique string of particle tree
    int m_parity;                                       ///< Intrinsic parity of particle
    int m_polState;                                     ///< polarisation state 
    unsigned int m_index;                               ///< Index, for constructing four-momenta
    unsigned int m_originalIndex;                       ///< Starting index, used in Bose-symmetrisation
    unsigned int m_orbital;                             ///< Orbital angular momentum between daughters
    unsigned int m_spinConfigurationNumber;             ///< Spin configuration quantum number
    unsigned int m_minL;                                ///< Minimum orbital angular momentum
    bool m_isHead;                                      ///< Flag that particle is head of decay chain
    bool m_isStateGood;                                 ///< Flag to check the decay is well-formed
    bool m_usesDefaultLineshape;                        ///< Flag to check if default shape is used
    std::vector<std::shared_ptr<Particle>> m_daughters; ///< Array of daughter particles
    std::vector<std::string> m_modifiers;               ///< Additional modifiers for amplitude

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
