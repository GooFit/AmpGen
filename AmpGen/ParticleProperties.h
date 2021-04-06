#ifndef AMPGEN_PARTICLEPROPERTIES_H
#define AMPGEN_PARTICLEPROPERTIES_H
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:04 GMT

#include <iostream>
#include <string>

#include "AmpGen/QuarkContent.h"
#include "AmpGen/Units.h"

namespace AmpGen
{
  /** @class ParticleProperties
    @brief Class that contains the PDG properties (mass, width, charges, etc.) for a single particle species, 

    Stores the information for a single particle species, the mass, the width, spin, various charges.
    The data is typically loaded from ParticlePropertiesList, which loads from the PDG provided 
    \code{cpp}
    mass_width.cvs
    \endcode  
    See that more details. 
    */


  class ParticleProperties
  {
    public:
      explicit ParticleProperties( const std::string& pdg_string = "" ); ///< Constructor from a string formatted by the PDG convention.  
      double mass()      const { return m_mass      * MeV; } ///< Returns mass of particle in MeV
      double mErrPlus()  const { return m_mErrPlus  * MeV; } ///< Returns +ve uncertainty on particle mass in MeV
      double mErrMinus() const { return m_mErrMinus * MeV; } ///< Returns -ve uncertainty on particle mass in MeV
      double width()     const { return m_width     * MeV; } ///< Returns width of particle in MeV
      double wErrPlus()  const { return m_wErrPlus  * MeV; } ///< Returns +ve uncertainty on particle width in MeV
      double wErrMinus() const { return m_wErrMinus * MeV; } ///< Returns -ve uncertainty on particle width in MeV
      double radius()    const;                              ///< Returns the effective interaction radius of the particle, i.e. for the Blatt-Weisskopf factors
      double lifetime()  const;                              ///< Returns the lifetime of the particle in ns

      int G()            const { return m_Gparity; }         ///< Returns the G-parity of the particle
      int P()            const { return m_parity;  }         ///< Returns the parity of the particle
      int C()            const { return m_Cparity; }         ///< Returns the C-parity of the particle
      int R()            const { return m_Rexist;  }         ///< Returns the R-parity of the particle
      int pdgID()        const { return m_pdgID;   }         ///< Returns the PDG id of the particle.
      int twoSpin()      const { return m_twoSpin; }         ///< Returns twice the spin of the particle 
      int charge()       const { return m_charge; }          ///< Returns the (electrical) charge of the particle
      char S()           const { return m_status; }          ///< Returns the existence status of the particle, i.e. whether it is confirmed by multiple experiments
      std::string I()    const { return m_isospin; }     ///< Returns the isospin of the particle as a string.
      std::string J()    const;     
      std::string label()    const { return m_texName; }     ///< Returns the LaTeX formatted label for the particle
      std::string name()     const;                          ///< Returns the particle name
      std::string spinName() const;                          ///< Returns the name of the particles spin. 
      bool isValid()         const { return m_isValid; }     ///< Check if the particle properties have been configured correctly 
      bool hasDistinctAnti() const;                          ///< Check if the particle has a distinct antiparticle
      bool isNonResonant()   const;                          ///< Check is this is a nonresonant `quasi-particle'
      bool isFermion()       const;                          ///< Check if the particle is a fermion, i.e. if the spin 1/2, 3/2, ...
      bool isBoson()         const;                          ///< Check if the particle is a boson, i.e. if the spin 0, 1, 2...
      bool isNeutrino()      const;                          ///< Check if the particle is a neutrino

      bool isPhoton()        const;                          ///< Check if the particle is a photon 
      void setProperty(const std::string& key, const std::string& value); ///< set a propery of a particle by key 
      const QuarkContent& quarkContent() const { return m_quarkContent; } ///< Returns the particle's quark content

      void setLabel( const std::string& label ) { m_texName = label; }          ///< Set the LaTeX label of the particle
      void setName( const std::string& name ) { 
        m_customName = true; 
        m_name = name; }                ///< Set the name of the particle

      void print( std::ostream& out = std::cout ) const;

      bool operator==( const ParticleProperties& rhs ) const;
      bool operator< ( const ParticleProperties& rhs ) const;
      bool operator> ( const ParticleProperties& rhs ) const;
      bool operator<=( const ParticleProperties& rhs ) const;
      bool operator>=( const ParticleProperties& rhs ) const;

      bool antiThis();                                                          ///< Change this particle to its antiparticle 
      ParticleProperties anti() const;                                          ///< Return the antiparticle of this particle

      static const ParticleProperties* get( const std::string& name, const bool& quiet=false ); 

    private:
      double m_mass{0};                  ///< mass [GeV]
      double m_mErrPlus{0};              ///< +ve mass error [GeV]
      double m_mErrMinus{0};             ///< -ve mass error [GeV]
      double m_width{0};                 ///< width [GeV]
      double m_wErrPlus{0};              ///< +ve width error [GeV]
      double m_wErrMinus{0};             ///< -ve width error [GeV]
      double m_radius{0};                ///< hadronic radius
      int m_Gparity{0};                  ///< G-parity
      int m_parity{0};                   ///< Parity
      int m_Cparity{0};                  ///< Charge 'parity'
      int m_pdgID{0};                    ///< PDG id
      int m_Rexist{0};                   ///< likelihood of existence, baryons only
      int m_charge{0};                   ///< electrical charge
      int m_twoSpin{0};                  ///< twice the spin 
      std::string m_isospin{""};         ///< isospin
      std::string m_name{""};            ///< particle name
      std::string m_texName{""};         ///< latex label of particle
      std::string m_chargeString{""};    ///< string for particle charge
      char m_Aformat;                    ///< anti-particle format character
      char m_status;                     ///< status (estalished or not etc)
      QuarkContent m_quarkContent;       ///< The quark content of the state (uD, uud etc.)
      bool m_isValid;                    ///< Flag to check whether the ParticleProperties have configured correctly 
      bool m_customName = {false};       ///< Flag to make custom name 
      void antiQuarks();
      void antiCharge();
      int chargeFromString( const std::string& ch, bool& status ) const;
  };
  std::ostream& operator<<( std::ostream& out, const AmpGen::ParticleProperties& pp );
}

#endif
//
