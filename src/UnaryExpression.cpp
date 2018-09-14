#include <algorithm>
#include <complex>
#include <memory>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpressionBase.h"

using namespace AmpGen;

DEFINE_CAST( Log )
DEFINE_CAST( Sqrt )
DEFINE_CAST( Exp )
DEFINE_CAST( Abs )
DEFINE_CAST( Sin )
DEFINE_CAST( Cos )
DEFINE_CAST( Tan )
DEFINE_CAST( ASin )
DEFINE_CAST( ACos )
DEFINE_CAST( ATan )
DEFINE_CAST( Norm )
DEFINE_CAST( Conj )
DEFINE_CAST( Real )
DEFINE_CAST( Imag )

std::complex<double> Log::operator()() const { return log( m_expression() ); }
std::complex<double> Sqrt::operator()() const { return sqrt( m_expression() ); }
std::complex<double> Exp::operator()() const { return exp( m_expression() ); }
std::complex<double> Abs::operator()() const { return std::abs( m_expression() ); }
std::complex<double> Sin::operator()() const { return sin( m_expression() ); }
std::complex<double> Cos::operator()() const { return cos( m_expression() ); }
std::complex<double> Tan::operator()() const { return tan( m_expression() ); }
std::complex<double> ASin::operator()() const { return asin( m_expression() ); }
std::complex<double> ACos::operator()() const { return acos( m_expression() ); }
std::complex<double> ATan::operator()() const { return atan( m_expression() ); }
std::complex<double> Norm::operator()() const { return std::norm( m_expression() ) ; }
std::complex<double> Conj::operator()() const { return std::conj( m_expression() ) ; }
std::complex<double> Real::operator()() const { return std::real( m_expression() ) ; }
std::complex<double> Imag::operator()() const { return std::imag( m_expression() ) ; }

Log::Log( const Expression& other ) : IUnaryExpression( other ) {}
Sqrt::Sqrt( const Expression& other ) : IUnaryExpression( other ) {}
Exp::Exp( const Expression& other ) : IUnaryExpression( other ) {}
Abs::Abs( const Expression& other ) : IUnaryExpression( other ) {}
Sin::Sin( const Expression& other ) : IUnaryExpression( other ) {}
Cos::Cos( const Expression& other ) : IUnaryExpression( other ) {}
Tan::Tan( const Expression& other ) : IUnaryExpression( other ) {}
ASin::ASin( const Expression& other ) : IUnaryExpression( other ) {}
ACos::ACos( const Expression& other ) : IUnaryExpression( other ) {}
ATan::ATan( const Expression& other ) : IUnaryExpression( other ) {}
Norm::Norm( const Expression& other ) : IUnaryExpression( other ) {}
Conj::Conj( const Expression& other ) : IUnaryExpression( other ) {}
Real::Real( const Expression& other ) : IUnaryExpression( other ) {}
Imag::Imag( const Expression& other ) : IUnaryExpression( other ) {}

std::string Log::to_string( const ASTResolver* resolver) const { return "log("       + m_expression.to_string(resolver) + ")"; }
std::string Sqrt::to_string(const ASTResolver* resolver) const { return "sqrt("      + m_expression.to_string(resolver) + ")"; }
std::string Exp::to_string( const ASTResolver* resolver) const { return "exp("       + m_expression.to_string(resolver) + ")"; }
std::string Abs::to_string( const ASTResolver* resolver) const { return "fabs("      + m_expression.to_string(resolver) + ")"; }
std::string Sin::to_string( const ASTResolver* resolver) const { return "sin("       + m_expression.to_string(resolver) + ")"; }
std::string Cos::to_string( const ASTResolver* resolver) const { return "cos("       + m_expression.to_string(resolver) + ")"; }
std::string Tan::to_string( const ASTResolver* resolver) const { return "tan("       + m_expression.to_string(resolver) + ")"; }
std::string ASin::to_string(const ASTResolver* resolver) const { return "asin("      + m_expression.to_string(resolver) + ")"; }
std::string ACos::to_string(const ASTResolver* resolver) const { return "acos("      + m_expression.to_string(resolver) + ")"; }
std::string ATan::to_string(const ASTResolver* resolver) const { return "atan("      + m_expression.to_string(resolver) + ")"; }
std::string Norm::to_string(const ASTResolver* resolver) const { return "std::norm(" + m_expression.to_string(resolver) + ")" ;}
std::string Conj::to_string(const ASTResolver* resolver) const { return "std::conj(" + m_expression.to_string(resolver) + ")" ;}
std::string Real::to_string(const ASTResolver* resolver) const { return "std::real(" + m_expression.to_string(resolver) + ")" ;}
std::string Imag::to_string(const ASTResolver* resolver) const { return "std::imag(" + m_expression.to_string(resolver) + ")" ;}

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
