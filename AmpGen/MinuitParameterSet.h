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
    MinuitParameterSet();
    MinuitParameterSet( const MinuitParameterSet& other );
    ~MinuitParameterSet();

    MinuitParameterSet getFloating();

    bool add( MinuitParameter* parPtr );
    MinuitParameter* add( const std::string& name, const unsigned int& flag, const double& mean, const double& sigma,
                          const double& min = 0, const double& max = 0 );
    bool unregister( MinuitParameter* patPtr );
    MinuitParameter* addOrGet( const std::string& name, const unsigned int& flag, const double& mean,
                               const double& sigma, const double& min = 0, const double& max = 0 );
    void loadFromStream();
    void loadFromFile( const std::string& name );
    void resetToInit();
    unsigned int size() const;

    MinuitParameter* getParPtr( unsigned int i ) const;

    std::map<std::string, MinuitParameter*>& map() { return _keyAccess; }
    const std::map<std::string, MinuitParameter*>& const_map() const { return _keyAccess; }
    std::vector<MinuitParameter*>::const_iterator cbegin() const { return _parPtrList.cbegin(); }
    std::vector<MinuitParameter*>::const_iterator cend()   const { return _parPtrList.cend(); }

    std::vector<MinuitParameter*> parPtrs() { return _parPtrList; }
    std::vector<MinuitParameter*>::iterator       begin() { return _parPtrList.begin(); }
    std::vector<MinuitParameter*>::iterator       end()   { return _parPtrList.end(); }
    std::vector<MinuitParameter*>::const_iterator begin() const { return _parPtrList.cbegin(); }
    std::vector<MinuitParameter*>::const_iterator end()   const { return _parPtrList.cend(); }

    void deleteListAndObjects();
    void deleteListKeepObjects();

    void print( std::ostream& os = std::cout ) const;
    void printVariable( std::ostream& os = std::cout ) const;

    void set( const MinuitParameterSet& mps );
    MinuitParameter* at( const std::string& key );
    MinuitParameter* operator[]( const std::string& key );
    MinuitParameter* operator[]( const std::string& key ) const;
    MinuitParameter* operator[]( const unsigned int& key );
    MinuitParameter* find( const std::string& key ) const 
    {
      auto it = _keyAccess.find(key);
      return it == _keyAccess.end() ? nullptr : it->second;   
    }
    double operator()( const std::string& name )
    {
      if ( _keyAccess.find( name ) == _keyAccess.end() ) {
        std::cout << "Cannot find parameter " << name << std::endl;
      }
      return _keyAccess[name]->mean();
    }
  private:
    void tryParameter( const std::vector<std::string>& line );
    void tryAlias( const std::vector<std::string>& line );
    std::vector<MinuitParameter*> _parPtrList;
    std::vector<MinuitExpression*> _aliasList;
    std::map<std::string, MinuitParameter*> _keyAccess;
    bool addToEnd( MinuitParameter* parPtr );

  };

} // namespace AmpGen
#endif
//
