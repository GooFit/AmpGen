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
  class ParticleProperties
  {
  private:
    double m_mass;              ///< mass [MeV]
    double m_mErrPlus;          ///< +ve mass error [MeV]
    double m_mErrMinus;         ///< -ve mass error [MeV]
    double m_width;             ///< width [MeV]
    double m_wErrPlus;          ///< +ve width error [MeV]
    double m_wErrMinus;         ///< -ve width error [MeV]
    double m_Radius;            ///< hadronic radius
    int m_Gparity;              ///< G-parity
    int m_Parity;               ///< Parity
    int m_Cparity;              ///< Charge 'parity'
    int m_pdgID;                ///< PDG id
    int m_Rexist;               ///< likelihood of existence, baryons only
    int m_charge;               ///< electrical charge
    int m_twoSpin;              ///< twice the spin 
    std::string m_Isospin;      ///< isospin
    std::string m_JtotalSpin;   ///< total spin
    std::string m_name;         ///< particle name
    std::string m_quarks;       ///< quark string
    std::string m_texName;      ///< latex label of particle
    std::string m_chargeString; ///< string for particle charge
    char m_Aformat;             ///< anti-particle format character
    char m_status;              ///< status (estalished or not etc)
    QuarkContent m_netQuarkContent;
    bool m_isValid;
    
    void setRadius();

    void antiQuarks();
    void antiQuarkContent();
    void antiCharge();
    int chargeFromString( const std::string& ch, bool& status ) const;

  public:
    double mass() const { return m_mass * MeV; }           ///< returns mass of particle in MeV
    double mErrPlus() const { return m_mErrPlus * MeV; }   ///< returns +ve uncertainty on particle mass in MeV
    double mErrMinus() const { return m_mErrMinus * MeV; } ///< returns -ve uncertainty on particle mass in MeV
    double width() const { return m_width * MeV; }         ///< returns width of particle in MeV
    double wErrPlus() const { return m_wErrPlus * MeV; }   ///< returns +ve uncertainty on particle width in MeV
    double wErrMinus() const { return m_wErrMinus * MeV; } ///< returns -ve uncertainty on particle width in MeV
    double radius() const;

    int G() const { return m_Gparity; }
    int P() const { return m_Parity; }
    int C() const { return m_Cparity; }
    int R() const { return m_Rexist; }
    int pdgID()    const { return m_pdgID; }
    int twoSpin()   const { return m_twoSpin ; }
    std::string I() const { return m_Isospin; }
    std::string J() const { return m_JtotalSpin; }
    int charge() const { return m_charge; }
    std::string label() const { return m_texName; }
    std::string quarks() const { return m_quarks; }
    std::string name() const;
    std::string chargeString() const { return m_chargeString; }
    std::string spinName() const; 

    char S() const { return m_status; }
    void setLabel( const std::string& label ) { m_texName = label; }
    void setName( const std::string& name ) { m_name = name; }
    const QuarkContent& netQuarkContent() const { return m_netQuarkContent; }
    
    bool isValid() const { return m_isValid; }

    bool hasDistinctAnti() const;
    bool barred() const;
    bool isItsOwnAnti() const { return !hasDistinctAnti(); }

    ParticleProperties( const std::string& pdg_string = "" );

    void print( std::ostream& out = std::cout ) const;

    bool operator==( const ParticleProperties& rhs ) const;
    bool operator< ( const ParticleProperties& rhs ) const;
    bool operator> ( const ParticleProperties& rhs ) const;
    bool operator<=( const ParticleProperties& rhs ) const;
    bool operator>=( const ParticleProperties& rhs ) const;

    bool antiThis();
    ParticleProperties anti() const;

    bool isNonResonant() const;
    bool isFermion() const ;
    bool isBoson() const;
    static const ParticleProperties* get( const std::string& name, const bool& quiet=false ); 
  };
} // namespace AmpGen
std::ostream& operator<<( std::ostream& out, const AmpGen::ParticleProperties& pp );

#endif
//
