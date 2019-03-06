#ifndef AMPGEN_PARTICLEPROPERTIESLIST_H
#define AMPGEN_PARTICLEPROPERTIESLIST_H
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:04 GMT

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/ParticleProperties.h"

namespace AmpGen
{
  class ParticlePropertiesList
  {
    static ParticlePropertiesList* ptr;
    std::map<int, std::pair<std::string, std::string>> m_latexLabels;
    std::string m_fname;
    std::vector<ParticleProperties> m_theList;
    std::map<std::string, ParticleProperties*> m_byName;
    std::map<int, ParticleProperties*> m_byID;
    double m_quasiStableThreshold;

    ParticlePropertiesList( const std::string& fname_in = "mass_width.csv" );

  protected:
    const std::vector<std::string> dirList() const;
    bool readFile( const std::string& fname );

  public:
    static const ParticlePropertiesList* getMe();
    static ParticlePropertiesList* getMutable();
    static const ParticleProperties* get( const std::string& name, const bool& quiet = false );
    static const ParticleProperties* get( const int& PDG, const bool& quiet = false );
    void makeAlias( const std::string& name, const std::string& alias );
    const ParticleProperties* find( const std::string& name, bool quiet = false ) const;
    const ParticleProperties* find( int pdg_id, bool quiet = false ) const;

    double quasiStableThreshold() const ;
    std::vector<std::string> getParticleNames() const;
    std::vector<int> getParticleIds() const;
    
    std::vector<ParticleProperties>::const_iterator begin() const { return m_theList.cbegin() ; }
    std::vector<ParticleProperties>::const_iterator   end() const { return m_theList.cend() ; }
    void print( std::ostream& out = std::cout ) const;
    bool readLatexLabels( const std::string& name );
    void makeMappings();
  };
  std::ostream& operator<<( std::ostream& out, const ParticlePropertiesList& ppl );
} // namespace AmpGen
#endif
