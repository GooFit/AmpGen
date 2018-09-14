#include "AmpGen/Spline.h"
#include "AmpGen/ASTResolver.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"

#include "TVectorD.h"

using namespace AmpGen; 

Expression AmpGen::Spline::operator()( const Expression& x , DebugSymbols* db )
{
  return Expression( SplineExpression(*this,x, db ) );
}

Expression AmpGen::getSpline( const std::string& name, const Expression& x, const std::string& arrayName,
    DebugSymbols* dbexpressions, const bool& continueSpline )
{
  double min, max, nBins( 0 );
  min   = AmpGen::NamedParameter<double>( name + "::Spline::Min", 0. ).getVal();
  max   = AmpGen::NamedParameter<double>( name + "::Spline::Max", 0. ).getVal();
  nBins = AmpGen::NamedParameter<double>( name + "::Spline::N"  , 0. ).getVal();
  std::string spline_name = name + "::Spline::"+arrayName;
  Spline shape( spline_name, nBins, min, max );
  return shape(x, dbexpressions); 
}

Expression Spline::eval(const Expression& xp, DebugSymbols* db )
{
  Expression x   = make_cse(xp);
  double spacing = ( m_max - m_min ) / ( (double)m_nKnots - 1. );
  Expression dx  = Fmod( x - m_min, spacing );
  Expression bin = ( x - m_min ) / spacing;
  Expression continuedValue            = 0;
  ADD_DEBUG( x   ,db );
  ADD_DEBUG( dx  ,db );
  ADD_DEBUG( bin ,db );
  ADD_DEBUG( m_points[bin], db );
  ADD_DEBUG( m_points[bin+1], db );
  ADD_DEBUG( m_points[bin+m_nKnots], db );
  ADD_DEBUG( m_points[bin+1+m_nKnots], db );

  Expression returnValue = Ternary( x > m_min && x < m_max,
      m_points[bin] + ( ( m_points[bin + 1] - m_points[bin] ) / spacing 
      - ( m_points[bin + 1 + m_nKnots] + 2 * m_points[bin+m_nKnots] ) * spacing / 6. ) * dx 
      + m_points[bin+m_nKnots] * dx * dx / 2. 
      + dx * dx * dx * ( m_points[bin+1+m_nKnots] - m_points[bin+m_nKnots] ) / ( 6. * spacing ),
      continuedValue );

  return SubTree(returnValue);
}
Spline::Spline(const std::string& name, const size_t& nBins, const double& min, const double& max, const std::vector<double>& values) :
  m_points( Parameter(name), 2*nBins ),
  m_name(name),
  m_nKnots(nBins),
  m_min(min),
  m_max(max),
  m_values(values){}

  void SplineTransfer::print() const { INFO( "Source: " << m_parameters[0]->name() ); }
  SplineTransfer::SplineTransfer() = default;

SplineTransfer::SplineTransfer( const SplineTransfer& other )
  : CacheTransfer()
  , m_transferMatrix( other.m_transferMatrix )
  , m_parameters( other.m_parameters )
  , m_min( other.m_min )
  , m_max( other.m_max )
    , m_address( other.m_address )
{
}

SplineTransfer::SplineTransfer( const unsigned int& address, const unsigned int& N, const double& min, const double& max )
  : CacheTransfer()
  , m_transferMatrix( TMatrixD( N - 2, N - 2 ) )
  , m_parameters( N, nullptr )
  , m_nKnots(N)
  , m_min( min )
  , m_max( max )
    , m_address( address )

{
  unsigned int size = N - 2;
  TMatrixD M( N - 2, N - 2 );
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
  for ( auto& x : m_parameters )
    if ( x == nullptr ) return false;
  return true;
}

void SplineTransfer::set( const unsigned int& N, AmpGen::MinuitParameter* f )
{
  m_parameters[N] = f;
}
void SplineTransfer::set( const unsigned int& N, const double& value )
{
  m_parameters[N] = new MinuitParameter("dumb",MinuitParameter::Fix,value,0);
}

void SplineTransfer::setAddress( const unsigned int& address ) { m_address = ( address ); }
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

void Spline::resolve( ASTResolver& resolver )
{
  m_points.m_address = resolver.addCacheFunction<SplineTransfer>(m_name,m_nKnots,m_min,m_max);
  auto splineTransfer = dynamic_cast<SplineTransfer*>( resolver.cacheFunctions[m_name].get() );
  if( resolver.mps == nullptr ) ERROR("Fix me!");
  if( m_values.size() == 0 ){
    for( unsigned int i = 0 ; i < m_nKnots; ++i ) 
      splineTransfer->set(i, resolver.mps->find(m_name+"::"+std::to_string(i)) );
  }
  else for( unsigned int i = 0 ; i < m_nKnots; ++i ) splineTransfer->set(i,m_values[i]) ;
}

SplineExpression::SplineExpression( const Spline& parent, const Expression& x , DebugSymbols* db) : 
  m_parent(std::make_shared<Spline>(parent)), 
  m_x(x), 
  m_eval( m_parent->eval(x, db) ) {}

void SplineExpression::resolve( ASTResolver& resolver ) 
{ 
  m_parent->resolve(resolver);  
  m_eval.resolve(resolver);
}

void Spline::set( const std::vector<double>& values )
{
  if( values.size() != m_nKnots ){ 
    ERROR("Sizes do not match");
    return;
  }
  m_values = values; 
}
