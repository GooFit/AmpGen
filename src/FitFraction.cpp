#include "AmpGen/FitFraction.h"

#include <complex>

#include "AmpGen/EventType.h"
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/FastIncoherentSum.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

FitFraction::FitFraction( const std::string& line, const AmpGen::EventType& evtType )
{
  auto tokens = split( line, ' ' );
  m_name      = tokens[1];
  m_value     = stod( tokens[2] );
  m_error     = stod( tokens[3] );
  if ( evtType.size() != 0 ) {
    std::vector<std::string> finalStates = evtType.finalStates();
  }
}

FFCalculator::FFCalculator( const std::string& name, 
                            FastCoherentSum* fcs, 
                            const std::vector<size_t>& indices,
                            const std::vector<size_t>& denom )
    : FFCalculator( name,fcs, indices, indices, denom ) {}

FFCalculator::FFCalculator( const std::string& name, 
                            FastCoherentSum* fcs, 
                            const std::vector<size_t>& indexI,
                            const std::vector<size_t>& indexJ, 
                            const std::vector<size_t>& denom )
    : index_i( indexI ), 
      index_j( indexJ ), 
      denom( denom ), 
      fcs( fcs ), 
      name( name ) {}

IFFCalculator::IFFCalculator( const size_t& index, FastIncoherentSum* fcs ) :
  index(index),
  fcs(fcs) {}

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

std::shared_ptr<Particle> FitFraction::particle() const { return std::make_shared<Particle>( m_name ); }

FitFraction::FitFraction( const std::string& name, const double& frac, const double& err )
    : m_name( name ), m_value( frac ), m_error( err )
{
}

