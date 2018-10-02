#include <algorithm>
#include <complex>
#include <memory>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpressionBase.h"

using namespace AmpGen;

DEFINE_UNARY_OPERATOR( Log , log )
DEFINE_UNARY_OPERATOR( Sqrt, sqrt )
DEFINE_UNARY_OPERATOR( Exp , exp )
DEFINE_UNARY_OPERATOR( Abs , abs )
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

Expression Log::d()  const { return 1. / arg(); }
Expression Sqrt::d() const { return 1. / ( 2 * fcn::sqrt( arg() ) ); }
Expression Exp::d()  const { return   fcn::exp(arg()) ; }
Expression Cos::d()  const { return - fcn::sin(arg()) ; } 
Expression Sin::d()  const { return   fcn::cos(arg()) ; }

Expression Tan::d()  const { return 1. / ( fcn::cos(arg()) * fcn::cos(arg()) );}

Expression ACos::d() const { return  1 / fcn::sqrt( 1 - arg()*arg() ); }
Expression ASin::d() const { return -1 / fcn::sqrt( 1 - arg()*arg() ) ; }
Expression ATan::d() const { return 1 / ( 1 + arg()*arg() ); }

Expression Conj::d() const { return 0;} 
Expression Abs::d()  const { return 0;} 
Expression Norm::d() const { return 0;} 
Expression Real::d() const { return 0;} 
Expression Imag::d() const { return 0;} 

void IUnaryExpression::resolve( ASTResolver& resolver ) { m_expression.resolve( resolver ); }
