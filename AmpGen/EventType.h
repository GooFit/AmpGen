#ifndef AMPGEN_EVENTTYPE_H
#define AMPGEN_EVENTTYPE_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

namespace AmpGen
{
  class Projection; 
  class Event; 
  /**@class EventType
   *  Deals with final state configuration of events,
   *  specifically dealing with the ordering of particles in trees.
   */

  class EventType
  {
  private:
    std::string m_mother;                            ///< name of decaying particle
    double m_motherMass;                             ///< mass of decaying particle
    std::vector<std::string> m_particleNames;        ///< names of decay products
    std::vector<std::string> m_particleNamesPickled; ///< names of decay product pickled for ROOT
    std::vector<double> m_particleMasses;            ///< masses of decay products
    bool m_timeDependent;

  public:
    EventType() = default;
    EventType( const std::vector<std::string>&, const bool& isTD = false );
    std::map<std::string, size_t> getEventFormat( const bool& outputNames = false ) const;

    std::pair<double, double> minmax( const std::vector<size_t>& indices, bool isGeV = false ) const;
    std::vector<std::vector<unsigned int>> getBosePairs() const;
    std::vector<double> masses() const;
    std::string mother() const;
    double mass( const size_t& index ) const;
    double motherMass() const;
    std::vector<std::string> finalStates() const;
    bool isTimeDependent() const;
    size_t eventSize() const;
    size_t size()      const;
    size_t dof()       const;
    std::string operator[]( const size_t& index ) const;

    std::string label( const size_t& index, bool isRoot = true ) const;
    std::string label( const std::vector<size_t>& index, bool isRoot = true ) const;
    std::vector<Projection> defaultProjections(const size_t& nBins) const;
    Projection projection(const size_t& nBins, const std::vector<size_t>& indices) const;

    bool operator==( const EventType& other ) const;
    bool has( const std::string& name ) const;
    EventType conj( const bool& headOnly = 0, const bool& dontConjHead = 0 ) const;
    
    void setMotherMass( const double& mass ){ m_motherMass = mass ; } 
    std::function<void( Event& )> symmetriser() const;
    std::pair<size_t, size_t> dim() const;
  };
  std::ostream& operator<<( std::ostream& os, const EventType& type );
} // namespace AmpGen

#endif
