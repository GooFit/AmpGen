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

MinuitParameterSet::MinuitParameterSet(const std::vector<MinuitParameter*>& params )
{
  for( auto& param : params ) add(param); 
}

MinuitParameterSet::MinuitParameterSet( const MinuitParameterSet& other )
  : m_parameters( other.m_parameters ), m_keyAccess( other.m_keyAccess ){}

MinuitParameterSet MinuitParameterSet::getFloating()
{
  MinuitParameterSet floating;
  for ( auto& param : *this ) {
    if ( param->isFree() && dynamic_cast<MinuitExpression*>(param) != nullptr ) 
      floating.add(param);
  }
  return floating;
}

bool MinuitParameterSet::addToEnd( MinuitParameter* parPtr )
{
  bool success = true;
  if ( nullptr == parPtr ) return false;
  m_parameters.push_back( parPtr );
  if ( m_keyAccess.find( parPtr->name() ) != m_keyAccess.end() ) {
    WARNING( "Parameter with name " << parPtr->name() << " already exists!" );
  }
  m_keyAccess[parPtr->name()] = parPtr;
  return success;
}

MinuitParameter* MinuitParameterSet::add( const std::string& name, const unsigned int& flag, const double& mean,
    const double& sigma, const double& min, const double& max )
{
  addToEnd( new MinuitParameter( name, Flag(flag), mean, sigma, min, max ) );
  return m_keyAccess[name];
}

bool MinuitParameterSet::add( MinuitParameter* parPtr ) { return addToEnd( parPtr ); }

bool MinuitParameterSet::unregister( MinuitParameter* parPtr )
{
  if ( m_parameters.end() == std::find( m_parameters.begin(), m_parameters.end(), parPtr ) ) {
    WARNING( "parPtr you want to unregister is not part of this list!" );
    return false;
  }
  m_parameters.erase( remove( m_parameters.begin(), m_parameters.end(), parPtr ), m_parameters.end() );
  return true;
}

unsigned int MinuitParameterSet::size() const { return m_parameters.size(); }

void MinuitParameterSet::print( std::ostream& os ) const
{
  for (size_t i = 0; i < size(); i++ ) {
    os << '\n' << *m_parameters[i];
  }
  os << '\n';
}
void MinuitParameterSet::printVariable( std::ostream& os ) const
{
  for (size_t i = 0; i < size(); i++ ) {
    if ( m_parameters[i]->flag() == Flag::Hide ) continue; /// hide parameter
    os << '\n' << *m_parameters[i];
  }
}

MinuitParameter* MinuitParameterSet::operator[]( const std::string& key )
{
  auto it = m_keyAccess.find( key );
  if ( it == m_keyAccess.end() ) {
    WARNING( "Parameter: " << key << " not found" );
  }
  return it->second;
}

MinuitParameter* MinuitParameterSet::operator[]( const std::string& key ) const
{
  auto it = m_keyAccess.find( key );
  if ( it == m_keyAccess.end() ) {
    WARNING( "Parameter: " << key << " not found" );
  }
  return it->second;
}

MinuitParameter* MinuitParameterSet::operator[]( const size_t& key ) { return m_parameters[key]; }

MinuitParameter* MinuitParameterSet::at( const std::string& key )
{
  if ( m_keyAccess.count( key ) == 0 ) {
    ERROR( key << " not found" );
    return nullptr;
  } else
    return m_keyAccess[key];
}

MinuitParameter* MinuitParameterSet::at( const size_t& index ) const
{
  if( index >= m_parameters.size() )
    ERROR( "Attempting to access parameter " << index << " when only " << m_parameters.size() << " have been defined" );
  return index < m_parameters.size() ? m_parameters[index] : nullptr; 
}


void MinuitParameterSet::tryParameter( const std::vector<std::string>& line )
{
  bool status = true;
  if ( line.size() == 4 || line.size() == 6 ) {
    bool hasLimits = line.size() == 6; 
    double mean = lexical_cast<double>( line[2], status );
    double step = lexical_cast<double>( line[3], status );
    double min  = hasLimits ? lexical_cast<double>( line[4], status ) : 0;
    double max  = hasLimits ? lexical_cast<double>( line[5], status ) : 0;
    if( !status ) return; 
    auto   flag = parse<Flag>( line[1] );
    if ( OptionsParser::printHelp() )
      INFO( "MINUIT: Registered " << line[0] << " (flag " << flag << ") = " << mean << ", step=" << step << " ("<< min << "," << max << ")" );
    add( new MinuitParameter( line[0], flag, mean, step, min, max ) ); 
  }
  if ( line.size() == 7 || line.size() == 11 ) {
    bool hasLimits = line.size() == 11;
    double mean_re = lexical_cast<double>( line[2], status );
    double step_re = lexical_cast<double>( line[3], status );
    double min_re  = hasLimits ? lexical_cast<double>( line[4], status ) : 0;
    double max_re  = hasLimits ? lexical_cast<double>( line[5], status ) : 0;
    double mean_im = lexical_cast<double>( line[5 + 2 *hasLimits], status );
    double step_im = lexical_cast<double>( line[6 + 2 *hasLimits], status );
    double min_im  = hasLimits ? lexical_cast<double>( line[9] , status ) : 0;
    double max_im  = hasLimits ? lexical_cast<double>( line[10], status ) : 0;

    if ( !status ) return;
    auto   flag_re = parse<Flag>(line[1]);
    auto   flag_im = parse<Flag>(line[4 + 2*hasLimits]);
    if ( OptionsParser::printHelp() ) {
      INFO( "MINUIT: Complex " << line[0] << "_Re (flag " << flag_re << ") = " << mean_re << ", step=" << step_re << " (" << min_re << "," << max_re << ")" );
      INFO( "MINUIT: Complex " << line[0] << "_Im (flag " << flag_im << ") = " << mean_im << ", step=" << step_im << " (" << min_im << "," << max_im << ")" );
    }
    add( new MinuitParameter( line[0] + "_Re", flag_re, mean_re, step_re, min_re, max_re ) );
    add( new MinuitParameter( line[0] + "_Im", flag_im, mean_im, step_im, min_im, max_im ) );
  }
}

void MinuitParameterSet::tryAlias( const std::vector<std::string>& line )
{
  if ( line.size() < 3 ) return;
  if ( line[1] == "=" ) {
    MinuitExpression* expr = new MinuitExpression(line, this );
    if ( expr->isGood() ) {
      addToEnd( expr );
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
}

void MinuitParameterSet::loadFromFile( const std::string& file )
{
  processFile( file, [this]( auto& line ) {
      this->tryParameter( split( line, {' ', '\t'} ) );
      this->tryAlias( split( line, {' ', '\t'} ) );
      } );
}

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

void MinuitParameterSet::rename(const std::string& name, const std::string& new_name)
{
  auto it = at(name);
  if( it == nullptr ){
    ERROR("Parameter: " << name << " not found");
    return;
  }
  if( name == new_name ) return;
  if( at(new_name) != nullptr ){
    ERROR("New key for " << name << " =  " << new_name << " already exists");
    return;
  } 
  it->setName(new_name);
  m_keyAccess.erase(name);
  m_keyAccess.emplace(new_name, it);
}

MinuitParameter* MinuitParameterSet::addOrGet( const std::string& name, const unsigned int& flag, const double& mean,
    const double& sigma, const double& min, const double& max )
{
  if ( m_keyAccess.count( name ) != 0 ) return m_keyAccess[name];
  return add( name, flag, mean, sigma, min, max );
}
MinuitParameterSet::const_iterator  MinuitParameterSet::cbegin() const { return m_parameters.cbegin(); }
MinuitParameterSet::const_iterator  MinuitParameterSet::cend()   const { return m_parameters.cend(); }
MinuitParameterSet::iterator        MinuitParameterSet::begin()        { return m_parameters.begin(); }
MinuitParameterSet::iterator        MinuitParameterSet::end()          { return m_parameters.end(); }
MinuitParameterSet::const_iterator  MinuitParameterSet::begin()  const { return m_parameters.cbegin(); }
MinuitParameterSet::const_iterator  MinuitParameterSet::end()    const { return m_parameters.cend(); }

MinuitParameter* MinuitParameterSet::find( const std::string& key ) const 
{
  auto it = m_keyAccess.find(key);
  return it == m_keyAccess.end() ? nullptr : it->second;   
}

double MinuitParameterSet::operator()( const std::string& name )
{
  if ( m_keyAccess.find(name) == m_keyAccess.end() ) {
    ERROR( "Cannot find parameter " << name );
  }
  return m_keyAccess[name]->mean();
}
