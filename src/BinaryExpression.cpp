#include <math.h>
#include <complex>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

DEFINE_BINARY_OPERATOR( Sum )
DEFINE_BINARY_OPERATOR( Sub )
DEFINE_BINARY_OPERATOR( Product )
DEFINE_BINARY_OPERATOR( Divide )
DEFINE_BINARY_OPERATOR( And )
DEFINE_BINARY_OPERATOR( GreaterThan )
DEFINE_BINARY_OPERATOR( LessThan )
DEFINE_BINARY_OPERATOR( Equal )
DEFINE_BINARY_OPERATOR( Pow )
DEFINE_BINARY_OPERATOR( Fmod )
DEFINE_BINARY_OPERATOR( ATan2 )

template < class condition > 
std::string bracketed( const Expression& expression, condition&& use_brackets, const ASTResolver* resolver=nullptr ){
  return use_brackets(expression) ? "(" + expression.to_string(resolver) +")" : expression.to_string(resolver);  
}

complex_t Sum::operator()()         const { return lval() + rval(); }
complex_t Sub::operator()()         const { return lval() - rval(); }
complex_t Product::operator()()     const { return lval() * rval(); }
complex_t Divide::operator()()      const { return lval() / rval(); }
complex_t And::operator()()         const { return complex_t(std::real(lval()) && std::real(rval()),0); }
complex_t GreaterThan::operator()() const { return complex_t(std::real(lval()) > std::real(rval()) ,0); }
complex_t LessThan::operator()()    const { return complex_t(std::real(lval()) < std::real(rval()) ,0); }
complex_t Pow::operator()()         const { return pow( lval(), rval() ); }
complex_t Fmod::operator()()        const { return 0; }
complex_t Equal::operator()()       const { return lval() == rval() ; } 
complex_t ATan2::operator()()       const { return atan2( std::real(lval() ), std::real(rval() ) ); }

std::string Sum::to_string(const ASTResolver* resolver)         const { 
  return lval.to_string(resolver) + " + " + rval.to_string(resolver) ;
}
std::string Sub::to_string(const ASTResolver* resolver)         const { 
  return lval.to_string(resolver) + "-"   + bracketed( rval, [](auto& expression){ return is<Sum>(expression) || is<Sub>(expression) ; } , resolver ) ; 
}
std::string Equal::to_string(const ASTResolver* resolver)       const { return "("     + lval.to_string(resolver) + " == "+ rval.to_string(resolver) +")"; }


std::string Product::to_string(const ASTResolver* resolver)     const { 
  auto use_brackets = [](auto& expression){ return is<Sum>(expression) || is<Sub>(expression); };
  return bracketed( lval,use_brackets,resolver) + "*" + bracketed(rval,use_brackets,resolver);
}

std::string Divide::to_string(const ASTResolver* resolver)      const { 
  //auto use_brackets = [](auto& expression){ return is<Sum>(expression) || is<Sub>(expression); };
  //auto use_brackets_r = [](auto& expression){ return is<IBinaryExpression>(expression) ; };
  return "(" + lval.to_string(resolver) + ")/("+ rval.to_string(resolver) +")";
  //return bracketed( lval, use_brackets,resolver) + "/"+ bracketed(rval, use_brackets_r,resolver);
}

std::string LessThan::to_string(const ASTResolver* resolver)    const { return "("     + lval.to_string(resolver) + "<"   + rval.to_string(resolver) +")"; }
std::string GreaterThan::to_string(const ASTResolver* resolver) const { return "("     + lval.to_string(resolver) + ">"   + rval.to_string(resolver) +")"; }
std::string And::to_string(const ASTResolver* resolver)         const { return "("     + lval.to_string(resolver) + "&&"  + rval.to_string(resolver) +")"; }
std::string Pow::to_string(const ASTResolver* resolver)         const { return "pow("  + lval.to_string(resolver) + ", "  + rval.to_string(resolver) +")"; }
std::string Fmod::to_string(const ASTResolver* resolver)        const { return "fmod(" + lval.to_string(resolver) + ","   + rval.to_string(resolver) +")"; }
std::string ATan2::to_string( const ASTResolver* resolver)      const { return "atan2("+ lval.to_string(resolver) + ","   + rval.to_string(resolver) +")"; }

void IBinaryExpression::resolve( ASTResolver& resolver ) const 
{
  lval.resolve( resolver );
  rval.resolve( resolver );
}

