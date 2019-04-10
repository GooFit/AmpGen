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
DEFINE_UNARY_OPERATOR( Abs , std::abs )
DEFINE_UNARY_OPERATOR( Sin , sin )
DEFINE_UNARY_OPERATOR( Cos , cos )
DEFINE_UNARY_OPERATOR( Tan , tan )
DEFINE_UNARY_OPERATOR( ASin, asin )
DEFINE_UNARY_OPERATOR( ACos, acos )
DEFINE_UNARY_OPERATOR( ATan, atan )
DEFINE_UNARY_OPERATOR( Norm, std::norm )
DEFINE_UNARY_OPERATOR( Conj, std::conj )
DEFINE_UNARY_OPERATOR( Real, std::real )
DEFINE_UNARY_OPERATOR( Imag, std::imag )
//DEFINE_UNARY_OPERATOR( LGamma, std::lgamma );
  //DEFINE_UNARY_OPERATOR( ISqrt, rsqrt )

ISqrt::ISqrt( const Expression& expression) : IUnaryExpression(expression) {} 
ISqrt::operator Expression() const { return Expression( std::make_shared<ISqrt>(*this) ) ; } 
complex_t ISqrt::operator()() const { return 1./sqrt( m_expression() ); } 
std::string ISqrt::to_string(const ASTResolver* resolver) const {   
  return resolver != nullptr && resolver->enableCuda()  ?
      "rsqrt("+m_expression.to_string(resolver)+")" :
    "1./sqrt("+m_expression.to_string(resolver)+")" ;
}

LGamma::LGamma( const Expression& expression) : IUnaryExpression(expression) {} 
LGamma::operator Expression() const { return Expression( std::make_shared<LGamma>(*this) ) ; }
complex_t LGamma::operator()() const { return std::lgamma( std::abs( m_expression() ) ); }
std::string LGamma::to_string(const ASTResolver* resolver) const {   
  return "std::lgamma(" + m_expression.to_string(resolver) + ")"; 
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
