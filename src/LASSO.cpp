#include <complex>
#include <string>
#include <vector>

#include "AmpGen/Factory.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

double LASSO::getVal() const
{
  double sum( 0 );
  for ( unsigned int i = 0; i < m_pdf->size(); ++i ) {
    std::complex<double> c_i = (*m_pdf)[i].coefficient;
    sum += sqrt( std::norm(c_i) * m_pdf->norm(i, i).real() );
  }
  return m_lambda * sum;
}

void LASSO::configure( const std::string& configString, 
                       const CoherentSum& pdf,
                       const MinuitParameterSet& mps )
{
  m_pdf       = &pdf;
  auto tokens = split( configString, ' ' );
  m_lambda    = stod( tokens[1] );
}

REGISTER( IExtendLikelihood, LASSO );
