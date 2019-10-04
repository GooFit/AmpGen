#include <memory.h>
#include <math.h>
#include <complex>
#include <iosfwd>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Utilities.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Simplify.h"
#include "AmpGen/ASTResolver.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

DEFINE_CAST(Constant )
DEFINE_CAST(Parameter )
DEFINE_CAST(SubTree )
DEFINE_CAST(Ternary )
DEFINE_CAST(Function )

Expression::Expression( const std::shared_ptr<IExpression>& expression ) : m_expression( expression ) {}

std::string  Expression::to_string(const ASTResolver* resolver) const
{
  if ( m_expression == nullptr || get() == nullptr ) {
    ERROR( "No expression contained in this node!" );
  }
  return m_expression->to_string(resolver);
}

IExpression* Expression::get() const { return m_expression.get(); }
complex_t    Expression::operator()() const { return ( *m_expression )(); }

template < class T > 
bool isEqual( const T& A, const double& B ){
  return std::abs(A-B) < std::numeric_limits<double>::epsilon();
}

bool isZero( const complex_t& A ){
  return std::abs(std::real(A)) < std::numeric_limits<double>::epsilon() && 
         std::abs(std::imag(A)) < std::numeric_limits<double>::epsilon();
}

std::string Constant::to_string(const ASTResolver* resolver) const {  
  auto rounded_string = [](const double& val ){
    std::string str = std::to_string (val);
    str.erase ( str.find_last_not_of('0') + 1, std::string::npos );  
    return str; 
  };
  std::string complex_type_string = resolver != nullptr && resolver->enableCuda() ? "ampgen_cuda::complex_t" : typeof<complex_t>() ;
  std::string literalSuffix = resolver != nullptr && resolver->enableCuda() ? "f" : ""; 
  return std::imag(m_value) == 0 ? "(" + rounded_string(std::real(m_value)) +literalSuffix + ")" : 
      complex_type_string +"("+rounded_string(std::real(m_value))+literalSuffix+","+rounded_string(std::imag(m_value))+literalSuffix+")";
}

Expression simplify_constant_addition( const Constant& constant, const Expression& expression )
{
  if( isZero( constant() ) )     return expression; 
  if( is<Constant>(expression) ) return Constant( constant() + expression() ); 
  if( is<Sum>(expression) ) return simplify_constant_addition( constant, cast<Sum>(expression).l() ) + cast<Sum>(expression).r();
  return Sum( constant, expression );
}

Expression simplify_product_addition( const Product& prod, const Expression& expression, bool& status )
{
  if( is<Constant>(prod.l()) && ( std::real( prod.l()() ) < 0.) ){
    return expression - Constant(-1.*complex_t(prod.l()()))*prod.r(); 
  }
  return Expression( Sum(expression, prod ) ); 
}

Expression AmpGen::operator+( const Expression& A, const Expression& B )
{
  bool status=true;
  if( is<Constant>(A) ) return simplify_constant_addition( cast<Constant>(A), B );
  if( is<Constant>(B) ) return simplify_constant_addition( cast<Constant>(B), A );
  if( is<Product>(A)  ){
    auto expr = simplify_product_addition( cast<Product>(A), B, status );
    if( status ) return expr;
  } 
  if( is<Product>(B)  ){
    return simplify_product_addition( cast<Product>(B), A , status );
  } 
  return Expression( Sum(A,B) );
}

Expression AmpGen::operator-( const Expression& A, const Expression& B )
{
  if ( is<Constant>(A) && isZero( A() ) )
    return -B;
  else if ( is<Constant>(B)&& isZero( B() ) )
    return A;
  else if ( is<Constant>(A) && is<Constant>(B) )
    return Constant( A() - B() );
  return Expression( Sub( A, B ) );
}

Expression simplify_multiplication( const Constant& constant, const Expression& expression )
{
  if( is<Constant>(expression ) ){   
    return Constant( constant() * expression() );
  }
  if( isEqual( constant(), 1 ) ) return expression;
  if( isZero( constant() ) ) return Constant(0);
  if( is<Product>(expression) )
  {
    return simplify_multiplication( constant , cast<Product>(expression).l() ) * cast<Product>(expression).r() ; /// always left associate multiplication 
  }
  return Product( constant, expression );
}

Expression simplify_sum_multiplication( const Sum& sum, const Expression& expression )
{
  return expression * sum.l() + expression * sum.r();
}

Expression simplify_sub_multiplication( const Sub& sum, const Expression& expression )
{
  return expression * sum.l() - expression * sum.r();
}

Expression AmpGen::operator*( const Expression& A, const Expression& B )
{
  if ( is<Constant>(A) ) return simplify_multiplication( cast<Constant>(A), B );
  if ( is<Constant>(B) ) return simplify_multiplication( cast<Constant>(B), A );
  if ( is<Divide>(A) )   return ( cast<Divide>(A).l() * B ) / ( cast<Divide>(A).r());
  if ( is<Divide>(B) )   return ( A * cast<Divide>(B).l() ) / ( cast<Divide>(B).r());
  if ( is<Sum>(A) || is<Sub>(A) ) return Expression( Product(B,A) );
  return Expression( Product( A, B ) );
}

Expression AmpGen::operator/( const Expression& A, const Expression& B )
{
  if( is<Constant>(A) && isEqual( A() , 0 ) )  return Constant(0);
  if( is<Constant>(B)){
    if(isEqual( B() , 1 ) ) return A;
    else return A * Constant( 1. / B() ) ; 
  }
  else if( is<Constant>(A) && is<Constant>(B) ) return Constant( A() / B() );
  else if( is<Product>(B)  ){
    auto as_prod = cast<Product>(B);
    if( is<Constant>( as_prod.l() ) ) return ( Constant( 1./as_prod.l()() ) * A )/ as_prod.r();
    if( is<Constant>( as_prod.r() ) ) return ( Constant( 1./as_prod.r()() ) * A )/ as_prod.l();
  }
  else if( is<Divide>(B) ) return ( A * cast<Divide>(B).r() ) / cast<Divide>(B).l();
  else if( is<Sqrt>(B) ) return ( A * fcn::isqrt( cast<Sqrt>(B).arg() ) ); 
  return Expression( Divide( A, B ) );
}
Expression AmpGen::operator&&( const Expression& A, const Expression& B ) { return Expression( And( A, B ) ); }

Expression AmpGen::operator==( const Expression& A, const Expression& B ){ return Equal(A,B) ; } 
Expression AmpGen::operator==( const double& A, const Expression& B ){ return Constant(A) == B ; } 
Expression AmpGen::operator==( const Expression& A, const double& B ){ return A == Constant(B) ; } 

Parameter::Parameter( const std::string& name, const double& defaultValue, const bool& resolved)
  : m_name( name )
  , m_defaultValue( defaultValue )
  , m_resolved( resolved )
{
}

std::string Parameter::to_string(const ASTResolver* resolver) const
{
  if ( m_resolved || resolver == nullptr ){
    return m_name; 
  }
  return resolver->resolvedParameter(this);
}

Expression AmpGen::operator<( const Expression& A, const Expression& B ) { return Expression( LessThan( A, B ) ); }
Expression AmpGen::operator>( const Expression& A, const Expression& B ) { return Expression( GreaterThan( A, B ) ); }

Expression Expression::operator-() const { return Constant( -1. ) * m_expression; }

Expression Expression::operator+=( const Expression& other ) const { return *this + other ; }
Expression Expression::operator-=( const Expression& other ) const { return *this - other ; }
Expression Expression::operator*=( const Expression& other ) const { return *this * other ; }
Expression Expression::operator/=( const Expression& other ) const { return *this /  other; }

std::ostream& AmpGen::operator<<( std::ostream& os, const Expression& expression ) {  return os << expression.to_string() ; } 
Expression::Expression( const double& value ) : m_expression( std::make_shared<Constant>( value ) ) {}
Expression::Expression( const complex_t& value ) : m_expression( std::make_shared<Constant>( value ) ) {}
Expression::Expression() : m_expression( std::make_shared<Constant>( 0. ) ) {}

void Expression::resolve( ASTResolver& resolver ) const { m_expression->resolve( resolver ); }
                           
void Constant::resolve( ASTResolver& resolver ) const {}

void Parameter::resolve( ASTResolver& resolver ) const
{
  if( !m_resolved ) resolver.resolve(*this);
}

Ternary::Ternary( const Expression& cond, const Expression& v1, const Expression& v2 )
  : m_cond( cond ), m_v1( v1 ), m_v2( v2 )
{
}
std::string Ternary::to_string(const ASTResolver* resolver) const
{
  return "(" + m_cond.to_string(resolver) + "?" + m_v1.to_string(resolver) + ":" + m_v2.to_string(resolver) + ")";
}

void Ternary::resolve( ASTResolver& resolver ) const
{
  m_cond.resolve( resolver );
  m_v1.resolve( resolver );
  m_v2.resolve( resolver );
}

SubTree::SubTree( const Expression& other ) : 
  m_expression( other ),
  m_key( FNV1a_hash(m_expression.to_string() ) ) {}

Function::Function( const std::string& name, const std::vector<Expression>& args ) :
  m_name(name),
  m_args(args)
{ 
}

std::string Function::to_string(const ASTResolver* resolver) const { 
  std::string rt = m_name + "(";
  if( m_args.size() == 0 ) return rt + ")";
  for( auto& arg : m_args ) rt += arg.to_string(resolver) + ", ";
  return rt.substr(0,rt.size()-2) + ")";
}
void Function::resolve( ASTResolver& resolver ) const {}

std::string SubTree::to_string(const ASTResolver* /*resolver*/) const 
{
  return "v"+ std::to_string(key());
}

void SubTree::resolve( ASTResolver& resolver ) const  
{ 
  resolver.resolve( *this );
  m_expression.resolve( resolver );
}

Expression AmpGen::make_cse( const Expression& A , bool simplify )
{
  if( is<Constant>(A) || is<Parameter>(A) || is<SubTree>(A) ) return A;
  SubTree cse = SubTree(A);
  if( !simplify ) return Expression(cse);
  auto ordered = Expression( NormalOrderedExpression(A) );
  cse.setKey( FNV1a_hash(ordered.to_string()) );
  return Expression(cse);
}

uint64_t SubTree::key() const { return m_key ; } // FNV1a_hash(m_expression.to_string()); }
void     SubTree::setKey( const size_t& new_key )
{
  m_key = new_key;
}
template <class TYPE> Expression simplify_constant_unary( const Expression& arg )
{
  return is<Constant>(arg) ? Expression(TYPE(arg)()) : Expression(TYPE(arg)); 
}

Expression AmpGen::fcn::sqrt(  const Expression& expression ) { return simplify_constant_unary<Sqrt>(expression) ; }
Expression AmpGen::fcn::abs(   const Expression& expression ) { return simplify_constant_unary<Abs>(expression); }
Expression AmpGen::fcn::sin(   const Expression& expression ) { return simplify_constant_unary<Sin>(expression); } 
Expression AmpGen::fcn::cos(   const Expression& expression ) { return simplify_constant_unary<Cos>(expression) ; } 
Expression AmpGen::fcn::tan(   const Expression& expression ) { return simplify_constant_unary<Tan>(expression) ; } 

Expression AmpGen::fcn::acos(   const Expression& expression ){ return simplify_constant_unary<ACos>(expression) ; } 
Expression AmpGen::fcn::asin( const Expression& expression){ return simplify_constant_unary<ASin>(expression) ; } 
Expression AmpGen::fcn::atan( const Expression& expression){ return simplify_constant_unary<ATan>(expression) ; } 
Expression AmpGen::fcn::atan2( const Expression& y, const Expression& x){ return ATan2(y,x) ; } 

Expression AmpGen::fcn::conj(  const Expression& expression ) { return is<Parameter>(expression) ? expression : simplify_constant_unary<Conj>(expression); } 
Expression AmpGen::fcn::norm(  const Expression& expression ) { return simplify_constant_unary<Norm>(expression); }
Expression AmpGen::fcn::exp(   const Expression& expression ) { return simplify_constant_unary<Exp>(expression) ; }
Expression AmpGen::fcn::log(   const Expression& expression ) { return simplify_constant_unary<Log>(expression) ; }

Expression AmpGen::fcn::safe_sqrt( const Expression& x )
{
  return Ternary( x > 0, fcn::sqrt(x), 0 );
}

Expression AmpGen::fcn::pow( const Expression& expression, const Expression& co )
{
  if( ! is<Constant>(co) ) return Pow(expression,co);
  if( isEqual( co(), 0 ) ) return 1;
  if( isEqual( co(), 1 ) ) return expression;
  else return Pow(expression, co ); 
}

Expression AmpGen::fcn::complex_sqrt( const Expression& expression )
{
  if( is<Constant>(expression ) ) return sqrt( expression() );
  auto st = make_cse(expression);
  return Ternary( st > 0, Sqrt(st), Constant(0,1)*Sqrt(-st) );
}

Expression AmpGen::fcn::isqrt( const Expression& expression )
{
  if( is<Constant>(expression ) ) return 1./sqrt( expression() );
  return ISqrt( expression );
}

Expression AmpGen::fcn::fpow( const Expression& x, const int& n){
  if( n < 0 ) return 1./fpow(x,-n);
  Expression rt = 1;
  for( int y=0;y<n;++y) rt = rt * x;
  return rt;
}
