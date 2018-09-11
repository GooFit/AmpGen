#include "AmpGen/FitFraction.h"

#include <algorithm>
#include <complex>

#include "AmpGen/EventType.h"
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/FastIncoherentSum.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

double FFCalculator::operator()()
{

  std::complex<double> J = 0;
  auto& amp              = *fcs;
  for ( auto& i : index_i ) {
    for ( auto& j : index_j ) {
      J += amp[i].coefficient * std::conj( amp[j].coefficient ) * fcs->norm( i, j );
    }
  }
  double norm = 0;
  if ( denom.size() != 0 ) {
    std::complex<double> N = 0;
    for ( auto& i : denom ) {
      for ( auto& j : denom ) {
        N += amp[i].coefficient * std::conj( amp[j].coefficient ) * fcs->norm( i, j );
      }
    }
    norm = std::real( N );
  } else
    norm = fcs->norm();
  return std::real( J ) / norm;
}

double IFFCalculator::operator()()
{
  fcs->transferParameters();
  double J = std::norm( ( *fcs )[index].coefficient ) * std::real( fcs->norm( index ) );
  return J / fcs->norm();
}

FitFraction::FitFraction( const std::string& line, const AmpGen::EventType& evtType )
{
  auto tokens = split( line, ' ' );
  m_name      = tokens[1];
  m_value     = stod( tokens[2] );
  m_error     = stod( tokens[3] );
  if ( evtType.size() != 0 ) {
    std::vector<std::string> finalStates = evtType.finalStates();
    // m_particle = std::make_shared<Particle>( tokens[1], finalStates );
  }
}

std::shared_ptr<Particle> FitFraction::particle() const { return std::make_shared<Particle>( m_name ); }

FitFraction::FitFraction( const std::string& name, const double& frac, const double& err )
    : m_name( name ), m_value( frac ), m_error( err )
{
}

FFCalculator::FFCalculator( const std::string& _name, FastCoherentSum* _fcs, const std::vector<unsigned int> indices,
                            const std::vector<unsigned int> _denom )
    : index_i( indices ), index_j( indices ), denom( _denom ), fcs( _fcs ), name( _name )
{
}

FFCalculator::FFCalculator( const std::string& _name, FastCoherentSum* _fcs, const std::vector<unsigned int> _indexI,
                            const std::vector<unsigned int> _indexJ, const std::vector<unsigned int> _denom )
    : index_i( _indexI ), index_j( _indexJ ), denom( _denom ), fcs( _fcs ), name( _name )
{
}
