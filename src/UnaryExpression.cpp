#include <memory.h>
#include <math.h>
#include <complex>
#include <string>
#include <utility>

#include "AmpGen/Expression.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Types.h"
#include "AmpGen/ASTResolver.h"

using namespace AmpGen;

template <class T>
T rsqrt( const T& arg ){ return 1. / sqrt(arg) ; } 

DEFINE_UNARY_OPERATOR( Log , log )
DEFINE_UNARY_OPERATOR( Sqrt, sqrt )
DEFINE_UNARY_OPERATOR( Exp , exp )
DEFINE_UNARY_OPERATOR( Sin , sin )
DEFINE_UNARY_OPERATOR( Cos , cos )
DEFINE_UNARY_OPERATOR( Tan , tan )
DEFINE_UNARY_OPERATOR( ASin, asin )
DEFINE_UNARY_OPERATOR( ACos, acos )
DEFINE_UNARY_OPERATOR( ATan, atan )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( Norm, std::norm )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( Real, std::real )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( Imag, std::imag )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( ISqrt, rsqrt )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( Conj, std::conj )
DEFINE_UNARY_OPERATOR_NO_RESOLVER( Abs , std::abs )

LGamma::LGamma( const Expression& expression) : IUnaryExpression(expression) {} 
LGamma::operator Expression() const { return Expression( std::make_shared<LGamma>(*this) ) ; }
complex_t LGamma::operator()() const { return std::lgamma( std::abs( m_expression() ) ); }
std::string LGamma::to_string(const ASTResolver* resolver) const {   
  return "std::lgamma(" + m_expression.to_string(resolver) + ")"; 
}

std::string ISqrt::to_string(const ASTResolver* resolver) const {   
  return resolver != nullptr && resolver->enableCuda()  ?
    "rsqrt("+m_expression.to_string(resolver)+")" :
    "1./sqrt("+m_expression.to_string(resolver)+")" ;
}

std::string Abs::to_string( const ASTResolver* resolver ) const
{
  return resolver != nullptr && resolver->enableAVX() ? 
    "abs(" + m_expression.to_string(resolver) +")" :
    "std::abs("+m_expression.to_string(resolver) +")";
}

std::string Conj::to_string( const ASTResolver* resolver ) const
{
  return resolver != nullptr && resolver->enableAVX() ? 
    "conj(" + m_expression.to_string(resolver) +")" :
    "std::conj("+m_expression.to_string(resolver) +")";
}

std::string Norm::to_string( const ASTResolver* resolver ) const
{
  return resolver != nullptr && resolver->enableAVX() ? 
    "norm(" + m_expression.to_string(resolver) +")" :
    "std::norm("+m_expression.to_string(resolver) +")";
}

std::string Real::to_string( const ASTResolver* resolver ) const
{
  return resolver != nullptr && resolver->enableAVX() ? 
    "real(" + m_expression.to_string(resolver) +")" :
    "std::real("+m_expression.to_string(resolver) +")";
}

std::string Imag::to_string( const ASTResolver* resolver ) const
{
  return resolver != nullptr && resolver->enableAVX() ? 
    "imag(" + m_expression.to_string(resolver) +")" :
    "std::imag("+m_expression.to_string(resolver) +")";
}


Expression Log::d()  const { return 1. / arg(); }
Expression Sqrt::d() const { return 1. / ( 2 * fcn::sqrt( arg() ) ); }
Expression Exp::d()  const { return   fcn::exp(arg()) ; }
Expression Cos::d()  const { return - fcn::sin(arg()) ; } 
Expression Sin::d()  const { return   fcn::cos(arg()) ; }

Expression Tan::d()  const { return 1. / ( fcn::cos(arg()) * fcn::cos(arg()) );}

Expression ACos::d() const { return  1 / fcn::sqrt( 1 - arg()*arg() ); }
Expression ASin::d() const { return -1 / fcn::sqrt( 1 - arg()*arg() ) ; }
Expression ATan::d() const { return 1 / ( 1 + arg()*arg() ); }

Expression ISqrt::d() const { return -1./ ( 2 * Pow( arg() , 1.5 ) ); }

Expression Conj::d() const { return 0;} 
Expression Abs::d()  const { return 0;} 
Expression Norm::d() const { return 0;} 
Expression Real::d() const { return 0;} 
Expression Imag::d() const { return 0;} 
Expression LGamma::d() const { return 0;}
void IUnaryExpression::resolve( ASTResolver& resolver ) const { m_expression.resolve( resolver ); }
