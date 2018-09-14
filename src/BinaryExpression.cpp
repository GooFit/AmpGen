#include <algorithm>
#include <complex>
#include <memory>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/ASTResolver.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

DEFINE_CAST( Sum )
DEFINE_CAST( Sub )
DEFINE_CAST( Product )
DEFINE_CAST( Divide )
DEFINE_CAST( And )
DEFINE_CAST( GreaterThan )
DEFINE_CAST( LessThan )
DEFINE_CAST( Equal )
DEFINE_CAST( Pow )
DEFINE_CAST( Fmod )

Sum::Sum( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
Sub::Sub( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
Product::Product( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
Divide::Divide( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
LessThan::LessThan( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
Pow::Pow( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
Fmod::Fmod( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
And::And( const Expression& l, const Expression& r ) : IBinaryExpression( l, r ) {}
GreaterThan::GreaterThan( const Expression& lval, const Expression& rval ) : IBinaryExpression( lval, rval ) {}
Equal::Equal( const Expression& lval, const Expression& rval ) : IBinaryExpression( lval, rval ) {}

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
  auto use_brackets = [](auto& expression){ return is<Sum>(expression) || is<Sub>(expression); };
  auto use_brackets_r = [](auto& expression){ return is<IBinaryExpression>(expression) ; };
  return bracketed( lval, use_brackets,resolver) + "/"+ bracketed(rval, use_brackets_r,resolver);
}

std::string LessThan::to_string(const ASTResolver* resolver)    const { return "("     + lval.to_string(resolver) + "<"   + rval.to_string(resolver) +")"; }
std::string GreaterThan::to_string(const ASTResolver* resolver) const { return "("     + lval.to_string(resolver) + ">"   + rval.to_string(resolver) +")"; }
std::string And::to_string(const ASTResolver* resolver)         const { return "("     + lval.to_string(resolver) + "&&"  + rval.to_string(resolver) +")"; }
std::string Pow::to_string(const ASTResolver* resolver)         const { return "pow("  + lval.to_string(resolver) + ", "  + rval.to_string(resolver) +")"; }
std::string Fmod::to_string(const ASTResolver* resolver)        const { return "fmod(" + lval.to_string(resolver) + ","   + rval.to_string(resolver) +")"; }

void IBinaryExpression::resolve( ASTResolver& resolver )
{
  lval.resolve( resolver );
  rval.resolve( resolver );
}

