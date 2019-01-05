// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/MinuitExpression.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

MinuitParameterSet::MinuitParameterSet() = default;

MinuitParameterSet::MinuitParameterSet( const MinuitParameterSet& other )
    : _parPtrList( other._parPtrList ), _keyAccess( other._keyAccess )
{
}

MinuitParameterSet MinuitParameterSet::getFloating()
{
  MinuitParameterSet floating;
  for ( auto& param : *this ) {
    if ( param->iFixInit() == MinuitParameter::Flag::Float ) floating.add(param);
  }
  return floating;
}

bool MinuitParameterSet::addToEnd( MinuitParameter* parPtr )
{
  bool success = true;
  if ( nullptr == parPtr ) return false;

  _parPtrList.push_back( parPtr );

  if ( _keyAccess.find( parPtr->name() ) != _keyAccess.end() ) {
    WARNING( "Parameter with name " << parPtr->name() << " already exists!" );
  }
  _keyAccess[parPtr->name()] = parPtr;
  return success;
}

MinuitParameter* MinuitParameterSet::add( const std::string& name, const unsigned int& flag, const double& mean,
                                          const double& sigma, const double& min, const double& max )
{
  addToEnd( new MinuitParameter( name, MinuitParameter::Flag(flag), mean, sigma, min, max ) );
  return _keyAccess[name];
}

bool MinuitParameterSet::add( MinuitParameter* parPtr ) { return addToEnd( parPtr ); }

bool MinuitParameterSet::unregister( MinuitParameter* parPtr )
{
  if ( _parPtrList.end() == std::find( _parPtrList.begin(), _parPtrList.end(), parPtr ) ) {
    WARNING( "parPtr you want to unregister is not part of this list!" );
    return false;
  }
  _parPtrList.erase( remove( _parPtrList.begin(), _parPtrList.end(), parPtr ), _parPtrList.end() );
  return true;
}

unsigned int MinuitParameterSet::size() const { return _parPtrList.size(); }

MinuitParameter* MinuitParameterSet::getParPtr( unsigned int i ) const
{
  if ( i >= _parPtrList.size() ) return nullptr;
  return _parPtrList[i];
}

void MinuitParameterSet::deleteListAndObjects()
{
  for ( std::vector<MinuitParameter*>::iterator it = _parPtrList.begin(); it != _parPtrList.end(); it++ ) {
    delete ( *it );
  }
  _parPtrList.clear();
}

void MinuitParameterSet::deleteListKeepObjects() { _parPtrList.clear(); }

void MinuitParameterSet::print( std::ostream& os ) const
{
  for ( unsigned int i = 0; i < size(); i++ ) {
    if ( nullptr == getParPtr( i ) ) continue;
    os << '\n';
    getParPtr( i )->print( os );
  }
  os << '\n';
}
void MinuitParameterSet::printVariable( std::ostream& os ) const
{
  for ( unsigned int i = 0; i < size(); i++ ) {
    if ( nullptr == getParPtr( i ) ) continue;
    if ( getParPtr( i )->iFixInit() == 1 ) continue; /// hide parameter
    os << '\n';
    getParPtr( i )->print( os );
  }
}

MinuitParameter* MinuitParameterSet::operator[]( const std::string& key )
{
  auto it = _keyAccess.find( key );
  if ( it == _keyAccess.end() ) {
    WARNING( "Parameter: " << key << " not found" );
  }
  return it->second;
}
MinuitParameter* MinuitParameterSet::operator[]( const std::string& key ) const
{
  auto it = _keyAccess.find( key );
  if ( it == _keyAccess.end() ) {
    WARNING( "Parameter: " << key << " not found" );
  }
  return it->second;
}
MinuitParameter* MinuitParameterSet::operator[]( const unsigned int& key ) { return _parPtrList[key]; }
MinuitParameter* MinuitParameterSet::at( const std::string& key )
{
  if ( _keyAccess.count( key ) == 0 ) {
    ERROR( key << " not found" );
    return nullptr;
  } else
    return _keyAccess[key];
}

void MinuitParameterSet::tryParameter( const std::vector<std::string>& line )
{
  if ( line.size() == 4 || line.size() == 6 ) {
    bool status = true;
    int flag    = lexical_cast<int>( line[1], status );
    double mean = lexical_cast<double>( line[2], status );
    double step = lexical_cast<double>( line[3], status );
    double min  = line.size() == 6 ? lexical_cast<double>( line[4], status ) : 0;
    double max  = line.size() == 6 ? lexical_cast<double>( line[5], status ) : 0;

    if ( status ) {
      if ( OptionsParser::printHelp() )
        INFO( "MINUIT: Registered " << line[0] << " (flag " << flag << ") = " << mean << ", step=" << step << " ("
                                    << min << "," << max << ")" );
      add( new MinuitParameter( line[0], MinuitParameter::Flag(flag), mean, step, min, max ) );
    }
    return;
  }
  if ( line.size() == 7  ) {

    bool status    = true;
    int flag_re    = lexical_cast<int>( line[1], status );
    double mean_re = lexical_cast<double>( line[2], status );
    double step_re = lexical_cast<double>( line[3], status );
    int flag_im    = lexical_cast<int>( line[4], status );
    double mean_im = lexical_cast<double>( line[5], status );
    double step_im = lexical_cast<double>( line[6], status );

    double min_re = 0;
    double max_re = 0;
    double min_im = 0;
    double max_im = 0;
    if ( !status ) return;

    if ( NamedParameter<unsigned int>( "MinuitParameterSet::RegulateParameters", 0 ) == 1 ) {
      // std::complex<double> z0 = mean_re*std::complex<double>( cos( mean_im) , sin(mean_im ) );
      if ( mean_re < 0 ) {
        mean_re = -mean_re;
        mean_im = mean_im + M_PI;
      }
      mean_im = atan2( sin( mean_im ), cos( mean_im ) );
      // std::complex<double> zP = mean_re*std::complex<double>( cos( mean_im) , sin(mean_im ) );
      // DEBUG("Normalised parameters = " << z0 << " " << zP );
      max_re = 1000;
    }
    if ( OptionsParser::printHelp() ) {
      INFO( "MINUIT: Complex " << line[0] << "_Re (flag " << flag_re << ") = " << mean_re << ", step=" << step_re
                               << " (" << min_re << "," << max_re << ")" );
      INFO( "MINUIT: Complex " << line[0] << "_Im (flag " << flag_im << ") = " << mean_im << ", step=" << step_im
                               << " (" << min_im << "," << max_im << ")" );
    }
    add( new MinuitParameter( line[0] + "_Re", MinuitParameter::Flag(flag_re), mean_re, step_re, min_re, max_re ) );
    add( new MinuitParameter( line[0] + "_Im", MinuitParameter::Flag(flag_im), mean_im, step_im, min_im, max_im ) );
  }
  if ( line.size() == 11  ) {

    bool status    = true;
    int flag_re    = lexical_cast<int>( line[1], status );
    double mean_re = lexical_cast<double>( line[2], status );
    double step_re = lexical_cast<double>( line[3], status );
    double min_re  = lexical_cast<double>( line[4], status );
    double max_re  = lexical_cast<double>( line[5], status );
    int flag_im    = lexical_cast<int>( line[6], status );
    double mean_im = lexical_cast<double>( line[7], status );
    double step_im = lexical_cast<double>( line[8], status );
    double min_im  = lexical_cast<double>( line[9], status );
    double max_im  = lexical_cast<double>( line[10], status );
    if ( !status ) return;

    add( new MinuitParameter( line[0] + "_Re", MinuitParameter::Flag(flag_re), mean_re, step_re, min_re, max_re ) );
    add( new MinuitParameter( line[0] + "_Im", MinuitParameter::Flag(flag_im), mean_im, step_im, min_im, max_im ) );
  }
}

void MinuitParameterSet::tryAlias( const std::vector<std::string>& line )
{
  if ( line.size() < 3 ) return;
  if ( line[1] == "=" ) {
    std::string name       = line[0];
    MinuitExpression* expr = new MinuitExpression( line, this );
    if ( expr->isGood() ) {
      _aliasList.push_back( expr );
      _keyAccess[name] = expr;
    } else {
      ERROR( "Expression is ill-formed: " << line[0] );
      delete expr;
    }
  }
}

void MinuitParameterSet::loadFromStream()
{
  auto ppfl = OptionsParser::getMe();
  
  std::vector<std::vector<std::string>> protoAliases;

  for ( auto it = ppfl->begin(); it != ppfl->end(); ++it ) {
    tryParameter( it->second );
    if ( it->second.size() >= 3 && it->second[1] == "=" ) protoAliases.push_back( it->second );
  }
  for ( auto& alias : protoAliases ) tryAlias( alias );
  //  print();
}

void MinuitParameterSet::loadFromFile( const std::string& file )
{
  processFile( file, [this]( auto& line ) {
    this->tryParameter( split( line, {' ', '\t'} ) );
    this->tryAlias( split( line, {' ', '\t'} ) );
  } );
}

MinuitParameterSet::~MinuitParameterSet() = default;

void MinuitParameterSet::set( const MinuitParameterSet& other )
{
  for ( auto& param : *this ) {
    auto otherValue = other[param->name()];
    if ( otherValue != nullptr ) param->setCurrentFitVal( otherValue->mean() );
  }
}

void MinuitParameterSet::resetToInit()
{
  for ( auto& param : *this ) param->resetToInit();
}

MinuitParameter* MinuitParameterSet::addOrGet( const std::string& name, const unsigned int& flag, const double& mean,
                                               const double& sigma, const double& min, const double& max )
{
  if ( _keyAccess.count( name ) != 0 ) return _keyAccess[name];
  return add( name, flag, mean, sigma, min, max );
}
