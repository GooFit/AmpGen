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
  class EventType; 
  std::ostream& operator<<( std::ostream& os, const EventType& type );

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
      std::map<std::string, unsigned> getEventFormat( const bool& outputNames = false ) const;

      /// Counts the number of particles in this event type with
      /// the same name as the index'th name.
      std::pair<unsigned,unsigned>  count(const unsigned& index) const;
      std::pair<double, double> minmax( const std::vector<unsigned>& indices) const;
      std::vector<double> masses() const;
      std::string mother() const;
      double mass( const unsigned& index ) const;
      double motherMass() const;
      std::vector<std::string> finalStates() const;
      bool isTimeDependent() const;
      unsigned eventSize() const;
      unsigned size()      const;
      unsigned dof()       const;
      std::string operator[]( const unsigned& index ) const;
      std::string decayDescriptor() const;
      std::string label( const unsigned& index, bool isRoot = true ) const;
      std::string label( const std::vector<unsigned>& index, bool isRoot = true ) const;
      std::vector<Projection> defaultProjections(const unsigned& nBins=100) const;
      Projection projection(const unsigned& nBins, const std::vector<unsigned>& indices, const std::string& observable = "mass2") const;

      bool operator==( const EventType& other ) const;
      bool has( const std::string& name ) const;

      void extendEventType( const std::string& branch );

      /// Functor to randomly symmetrise data of this event type, using the Fisher-Yates shuffle.
      std::function<void( Event& )> symmetriser() const;
      std::function<bool( Event&, const std::vector<int>& ids)> automaticOrdering() const;

      /// Calculates the number of spin indices associated with the initial and final state, i.e. the rank of the relevant transition matrix.
      std::pair<unsigned, unsigned> dim() const;
  
      friend std::ostream& AmpGen::operator<<( std::ostream& os, const EventType& type );

    private:
      std::string               m_mother;               ///< name of decaying particle
      double                    m_motherMass;           ///< mass of decaying particle
      std::vector<std::string>  m_particleNames;        ///< names of decay products
      std::vector<std::string>  m_particleNamesPickled; ///< names of decay product pickled for ROOT
      std::vector<double>       m_particleMasses;       ///< masses of decay products
      std::vector<bool>         m_ignore;               ///< Flag to ignore particle when reading events, for invisible or tags 
      bool                      m_timeDependent;        ///< Flag to include a decay time as the last element in the event vector
      std::vector<std::string>  m_eventTypeExtensions;  ///< extended event data
      std::pair<unsigned, unsigned> m_dim;              ///< Rank of the relevant transition matrix
      bool                      m_alt_part_names;       ///< alternative naming in ouput tree (e.g. Xi- pi+ pi+ becomes Xim pip0 pip1 rather than _1_Xi# _2_pi~ _3_pi~)
  };
} // namespace AmpGen

#endif
