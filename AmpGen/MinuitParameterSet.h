#ifndef AMPGEN_MINUITPARAMETERSET_H
#define AMPGEN_MINUITPARAMETERSET_H

// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT

#include <iostream>
#include <map>
#include <vector>

#include "AmpGen/MinuitParameter.h"

namespace AmpGen
{
  class MinuitExpression;

  class MinuitParameterSet
  {
  public:
    typedef std::vector<MinuitParameter*>::iterator iterator; 
    typedef std::vector<MinuitParameter*>::const_iterator const_iterator; 
    
    MinuitParameterSet();
    MinuitParameterSet(const std::vector<MinuitParameter*>& params );
    MinuitParameterSet( const MinuitParameterSet& other );
    ~MinuitParameterSet() = default;

    MinuitParameterSet getFloating();

    bool add( MinuitParameter* parPtr );
    MinuitParameter* add(const std::string& name, const Flag& flag, const double& mean, const double& sigma, const double& min = 0, const double& max = 0 );
    bool unregister( MinuitParameter* patPtr );
    MinuitParameter* addOrGet(const std::string& name, const Flag& flag, const double& mean,
                              const double& sigma, const double& min = 0, const double& max = 0 );
    void loadFromStream();
    void loadFromFile( const std::string& name );
    void resetToInit();
    void print( std::ostream& os = std::cout ) const;
    void printVariable( std::ostream& os = std::cout ) const;
    void set( const MinuitParameterSet& mps );
    void rename(const std::string& name, const std::string& new_name); 
    unsigned int size() const;

    const_iterator cbegin() const;
    const_iterator cend()   const;
    iterator       begin();
    iterator       end();
    const_iterator begin() const;
    const_iterator end()   const;
    
    MinuitParameter* at( const std::string& key );
    MinuitParameter* at( const size_t& index ) const;
    MinuitParameter* operator[]( const std::string& key );
    MinuitParameter* operator[]( const std::string& key ) const;
    MinuitParameter* operator[]( const size_t& key );
    MinuitParameter* find( const std::string& key ) const;
    double operator()( const std::string& name );
  private:
    void tryParameter( const std::vector<std::string>& line );
    void tryAlias( const std::vector<std::string>& line );
    bool addToEnd( MinuitParameter* parPtr );
 
    std::vector<MinuitParameter*>           m_parameters;
    std::map<std::string, MinuitParameter*> m_keyAccess;
  };
} // namespace AmpGen
#endif
//
