#include "AmpGen/FitFraction.h"

#include <complex>

#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

FitFraction::FitFraction( const std::string& line )
{
  auto tokens = split( line, ' ' );
  m_name      = tokens[1];
  m_value     = stod( tokens[2] );
  m_error     = stod( tokens[3] );
}

std::shared_ptr<Particle> FitFraction::particle() const { return std::make_shared<Particle>( m_name ); }

FitFraction::FitFraction( const std::string& name, const double& frac, const double& err )
  : m_name( name ), m_value( frac ), m_error( err )
{
}
void FitFraction::setFracErr(const double& f, const double& e)
{
  m_value = f;
  m_error = e;
}
double FitFraction::val() const { return m_value; }
double FitFraction::err() const { return m_error; }
std::string FitFraction::name() const { return m_name; }

bool AmpGen::operator  <(const FitFraction& lhs, const FitFraction& rhs){
  return std::abs(lhs.val()) < std::abs(rhs.val());
}
bool AmpGen::operator  >(const FitFraction& lhs, const FitFraction& rhs){
  return std::abs(lhs.val()) > std::abs(rhs.val());
}
bool AmpGen::operator ==(const FitFraction& lhs, const FitFraction& rhs){
  return lhs.name() == rhs.name();
}

std::ostream& AmpGen::operator <<(std::ostream& os, const FitFraction& obj )
{
  return os << std::left  << std::setw(60) << obj.name()              << " = " 
            << std::right << std::setw(8)  << round(obj.val()*100, 3) << " ± " 
            << std::left  << std::setw(8)  << round(obj.err()*100, 3) << " %";
}
