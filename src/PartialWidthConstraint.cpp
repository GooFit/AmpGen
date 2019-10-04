#include <complex>
#include <ostream>
#include <string>
#include <vector>

#include "AmpGen/Factory.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"
#include "AmpGen/MinuitParameterSet.h"

using namespace AmpGen;

double PartialWidthConstraint::getVal() const
{

  std::complex<double> numerator( 0, 0 );
  std::complex<double> denom( 0, 0 );
  for ( auto& ip : m_numComponents ) {
    for ( auto& jp : m_numComponents ) {
      numerator += ( *m_pdf )[ip].coefficient * std::conj( ( *m_pdf )[jp].coefficient ) * m_pdf->norm( ip, jp );
    }
  }

  for ( auto& ip : m_denComponents ) {
    for ( auto& jp : m_denComponents ) {
      denom += ( *m_pdf )[ip].coefficient * std::conj( ( *m_pdf )[jp].coefficient ) * m_pdf->norm( ip, jp );
    }
  }
  double ratio = std::real( numerator ) / std::real( denom );
  DEBUG( "Ratio of BRs = " << ratio << " weight = " << m_weight << " ratio to aim for = " << m_ratio );
  return m_weight * ( ratio - m_ratio ) * ( ratio - m_ratio );
}

void PartialWidthConstraint::configure( const std::string& configString, 
                                        const CoherentSum& pdf,
                                        const MinuitParameterSet& mps )
{

  m_pdf                                = &( pdf );
  auto tokens                          = split( configString, ' ' );
  const std::string name               = tokens[1];
  m_weight                             = stod( tokens[2] );
  m_ratio                              = stod( tokens[3] );
  std::vector<std::string> denChannels = NamedParameter<std::string>( name + "_denChannels" ).getVector();
  std::vector<std::string> numChannels = NamedParameter<std::string>( name + "_numChannels" ).getVector();
  for ( auto& p : denChannels ) m_denComponents.push_back( findIndex(pdf.matrixElements(), p ) );
  for ( auto& p : numChannels ) m_numComponents.push_back( findIndex(pdf.matrixElements(), p ) );
  
  INFO( "Constraining ratio of " );
  for ( unsigned int i = 0; i < m_denComponents.size(); ++i )
    INFO( denChannels[i] << " index = " << m_denComponents[i] );
  INFO( " to " );
  for ( unsigned int i = 0; i < m_numComponents.size(); ++i )
    INFO( numChannels[i] << " index = " << m_numComponents[i] );
}

REGISTER( IExtendLikelihood, PartialWidthConstraint );
