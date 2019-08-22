#ifndef AMPGEN_EVENTTYPE_H
#define AMPGEN_EVENTTYPE_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>
#include <initializer_list>

namespace AmpGen
{
  class Projection; 
  class Event; 
  /**@class EventType
     Deals with final state configuration of events,
     specifically dealing with the ordering of particles in trees.
   */

  class EventType
  {
    public:
      /// Default constructor
      EventType() = default;

      /// Takes a list of particles, beginning with the head of the decay and 
      /// then then final state particles, and a flag as to whether to include time dependence. 
      EventType( const std::vector<std::string>&, const bool& isTD = false );
     
      /// Returns the CP-conjugated event type. Can also require that only the initial/
      /// final state is conjugated. By default, conjugates both.   
      EventType conj( const bool& headOnly = 0, const bool& dontConjHead = 0 ) const;
      
      /// Returns the event format, that matches between expressions and the names used in Particle.  
      std::map<std::string, size_t> getEventFormat( const bool& outputNames = false ) const;

      /// Counts the number of particles in this event type with 
      /// the same name as the index'th name. 
      std::pair<size_t,size_t>  count(const size_t& index) const;
      std::pair<double, double> minmax( const std::vector<size_t>& indices, bool isGeV = false ) const;
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
      std::string decayDescriptor() const; 
      std::string label( const size_t& index, bool isRoot = true ) const;
      std::string label( const std::vector<size_t>& index, bool isRoot = true ) const;
      std::vector<Projection> defaultProjections(const size_t& nBins) const;
      Projection projection(const size_t& nBins, const std::vector<size_t>& indices, const std::string& observable = "mass2") const;

      bool operator==( const EventType& other ) const;
      bool has( const std::string& name ) const;

      void extendEventType( const std::string& branch );

      /// Functor to randomly symmetrise data of this event type, using the Fisher-Yates shuffle.  
      std::function<void( Event& )> symmetriser() const;

      /// Calculates the number of spin indices associated with the initial and final state, i.e. the rank of the relevant transition matrix. 
      std::pair<size_t, size_t> dim() const;

    private:
      std::string               m_mother;               ///< name of decaying particle
      double                    m_motherMass;           ///< mass of decaying particle
      std::vector<std::string>  m_particleNames;        ///< names of decay products
      std::vector<std::string>  m_particleNamesPickled; ///< names of decay product pickled for ROOT
      std::vector<double>       m_particleMasses;       ///< masses of decay products
      bool                      m_timeDependent;        ///< Flag to include a decay time as the last element in the event vector
      std::vector<std::string>  m_eventTypeExtensions;  ///< extended event data
      std::pair<size_t, size_t> m_dim;                  ///< Rank of the relevant transition matrix 
  };
  std::ostream& operator<<( std::ostream& os, const EventType& type );
} // namespace AmpGen

#endif
