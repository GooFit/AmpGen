#include "AmpGen/Spline.h"

#include <ostream>

#include "AmpGen/ASTResolver.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/MsgService.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include "TMatrixTUtils.h"
#include "TVectorT.h"

using namespace AmpGen; 

Spline::Spline(const std::string& name, 
    const size_t& nKnots, 
    const double& min, 
    const double& max ) :
  m_points( Parameter(name), 2*nKnots ),
  m_name(name),
  m_nKnots(nKnots),
  m_min(min),
  m_max(max) {}

Spline::Spline(const Spline& spline, 
    const Expression& x )
  : 
    m_points( spline.m_points ),
    m_name(   spline.m_name),
    m_nKnots( spline.m_nKnots),
    m_min(    spline.m_min ),
    m_max(    spline.m_max ),
    m_x  ( x ),
    m_eval( eval() )  {
    }
Expression Spline::operator()( const Expression& x )
{
  return Spline(*this,x); 
}

Expression AmpGen::getSpline( const std::string& name, const Expression& x, const std::string& arrayName,
    DebugSymbols* dbexpressions, const bool& continueSpline )
{
  double min, max, nBins( 0 );
  auto spline_params = NamedParameter<double>( name + "::Spline").getVector();
  if( spline_params.size() == 3 ){
    nBins = size_t( spline_params[0] );
    min   =         spline_params[1] ; 
    max   =         spline_params[2];
  }
  else {
    nBins = NamedParameter<double>( name + "::Spline::N"  , 0. );
    min   = NamedParameter<double>( name + "::Spline::Min", 0. );
    max   = NamedParameter<double>( name + "::Spline::Max", 0. );
  }
  std::string spline_name = name + "::Spline::"+arrayName;
  return Spline( spline_name, nBins, min,max )(x);
}

Expression Spline::eval() const 
{
  Expression x   = make_cse(m_x);
  double spacing = ( m_max - m_min ) / ( (double)m_nKnots - 1. );
  Expression dx  = Fmod( x - m_min, spacing );
  Expression bin = ( x - m_min ) / spacing;
  Expression continuedValue            = 0;

  Expression returnValue = Ternary( x > m_min && x < m_max,
      m_points[bin] + ( ( m_points[bin + 1] - m_points[bin] ) / spacing 
        - ( m_points[bin + 1 + m_nKnots] + 2 * m_points[bin+m_nKnots] ) * spacing / 6. ) * dx 
      + m_points[bin+m_nKnots] * dx * dx / 2. 
      + dx * dx * dx * ( m_points[bin+1+m_nKnots] - m_points[bin+m_nKnots] ) / ( 6. * spacing ),
      continuedValue );
  return make_cse(returnValue);
}

void SplineTransfer::print() const { INFO( "Source: " << m_parameters[0]->name() ); }
SplineTransfer::SplineTransfer() = default;

SplineTransfer::SplineTransfer( const SplineTransfer& other )
  : CacheTransfer(other.m_address, other.m_value, other.m_size)
  , m_transferMatrix( other.m_transferMatrix )
  , m_parameters( other.m_parameters )
  , m_min( other.m_min )
  , m_max( other.m_max )
{
}

SplineTransfer::SplineTransfer( const size_t& address, const unsigned int& N, const double& min, const double& max )
  : CacheTransfer(address)
  , m_transferMatrix( TMatrixD( N - 2, N - 2 ) )
  , m_parameters( N, nullptr )
  , m_nKnots(N)
  , m_min( min )
  , m_max( max )
{
  unsigned int size = N - 2;
  TMatrixD M(size, size);
  for ( unsigned int i = 0; i < size; ++i ) {
    M[i][i] = 4;
    if ( i != size - 1 ) {
      M[i][i + 1] = 1;
      M[i + 1][i] = 1;
    }
  }
  DEBUG( m_transferMatrix.GetNrows() << " x " << m_transferMatrix.GetNcols() << "  ; " << M.GetNrows() << " x "
      << M.GetNcols() );

  m_transferMatrix = M.Invert();
  DEBUG( "Leaving constructor ..." );
}

bool SplineTransfer::isConfigured()
{
  return std::all_of( m_parameters.begin(), m_parameters.end(), [](auto& p ){ return p != nullptr ; } );
}

void SplineTransfer::set( const unsigned int& N, MinuitParameter* f )
{
  m_parameters[N] = f;
}
void SplineTransfer::set( const unsigned int& N, const double& value )
{
  m_parameters[N] = new MinuitParameter("dumb", Flag::Fix, value, 0);
}

void SplineTransfer::transfer( CompiledExpressionBase* destination )
{
  unsigned int size = m_parameters.size() - 2;
  TVectorD L( size );
  for ( unsigned int i = 0; i < size; ++i ) {
    L[i] = m_parameters[i + 2]->mean() - 2 * m_parameters[i + 1]->mean() + m_parameters[i]->mean();
  }
  DEBUG( L.GetNrows() << "   " << m_transferMatrix.GetNrows() << "   " << m_transferMatrix.GetNcols() );
  auto mtv = m_transferMatrix * L;
  std::vector<double> mvectors( m_parameters.size(), 0 );
  double spacing = ( m_max - m_min ) / double( m_parameters.size() - 1);

  for ( int i = 0; i < mtv.GetNrows(); ++i ) mvectors[i + 1] = 6 * mtv[i] / ( spacing * spacing );

  for ( unsigned int i = 0; i < m_parameters.size(); ++i ) {
    DEBUG( "Knot["<<i<<"] : value = (" << m_parameters[i]->mean() << ", "  << m_address +i << ") "
        << " curvature = (" << mvectors[i] << " , " << m_address + m_nKnots +i << ")" );
    destination->setExternal( m_parameters[i]->mean(), m_address + i );
    destination->setExternal( mvectors[i], m_address + m_nKnots + i );  
  }
}

void Spline::resolve( ASTResolver& resolver ) const
{
  resolver.resolve(*this);
  m_x.resolve(resolver);
  m_eval.resolve(resolver);
}

std::string Spline::to_string(const ASTResolver* resolver) const { return m_eval.to_string(resolver) ; } 
Spline::operator Expression() { return Expression( std::make_shared<Spline>( *this ) ); }
complex_t Spline::operator()() const  { return 0; }
