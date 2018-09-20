#ifndef AMPGEN_EXPRESSION_H
#define AMPGEN_EXPRESSION_H

/// \defgroup ExpressionEngine Expressions 
/// \brief 
/// are a set of classes for constructing functions from primatives 
/// such as the binary operators \f$(+,-,/,\times)\f$, and elementary functions such as trigonometric 
/// functions and sqrts. 
///
/// AmpGen Expressions are a set of classes for constructing functions from primatives 
/// such as the binary operators \f$(+,-,/,\times)\f$, and elementary functions such as trigonometric 
/// functions and sqrts. 
/// These are then compiled to conventional C++ at runtime, to allow fast and flexible 
/// evalulation of amplitudes. 
/// Expression trees can be built using overloaded operators, so 
///
/// \code{.cpp}
/// Expression x = a+b; 
/// \endcode
///
/// will given an expression that can calculate the sum of expressions a and b 
/// [n.b. the auto keyword should be avoided if dealing with unary functions explicitly, but not if handled via the fcn:: namespace]
/// \author T.Evans


/// \ingroup ExpressionEngine macro ADD_DEBUG
/// Make a (named) debugging expression and add to a set of DebugSymbols.
#define ADD_DEBUG( X, Y )                                               \
  if ( Y != 0 ) Y->push_back( DebugSymbol( std::string( #X ), X ) );    \

/// \ingroup ExpressionEngine macro DEFINE_CAST 
/// Define a cast from an expression implementation to the shared_ptr wrapper. 
#define DEFINE_CAST( X )                                    \
  X::operator Expression() const                            \
{ return Expression( std::make_shared<X>( *this ) ); }  

/// \ingroup ExpressionEngine macro DECLARE_UNARY_OPERATOR
/// Macro to declare a unary operator, \ref ExpressionEngine "see IUnaryExpression"
#define DECLARE_UNARY_OPERATOR( X )                         \
  class X : public IUnaryExpression {                       \
    public:                                                 \
    X( const Expression& other );                           \
    virtual std::string to_string(const ASTResolver* resolver=nullptr) const override;         \
    virtual Expression d() const override;                  \
    operator Expression() const;                            \
    Expression clone() const override { return X(m_expression.clone());} \
    virtual std::complex<double> operator()() const override;          \
  }

/// \ingroup ExpressionEngine macro DECLARE_BINARY_OPERATOR
/// Macro to declare a binary operator, \ref ExpressionEngine "see IBinaryExpression"
#define DECLARE_BINARY_OPERATOR( X )                        \
  class X : public IBinaryExpression {                      \
    public:                                                 \
    X( const Expression& l, const Expression& r );          \
    virtual std::string to_string(const ASTResolver* resolver=nullptr) const override ;        \
    operator Expression() const ;                           \
    Expression clone() const override { return X(lval.clone(),rval.clone());} \
    virtual std::complex<double> operator()() const override;          \
  }

#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <functional>
#include "AmpGen/MsgService.h"
#include "AmpGen/MetaUtils.h"
namespace AmpGen
{
  struct ASTResolver;
  struct Expression;
  struct Parameter;

  typedef std::pair<std::string, Expression> DebugSymbol;
  typedef std::vector<DebugSymbol> DebugSymbols;

  /// \ingroup ExpressionEngine class IExpression 
  ///  Virtual base class for other expression tree components.
  ///  Implementations must permit the following operations on nodes.
  
  class IExpression {
    public:
      /// Called to convert the Expression tree into source code.
      /// \return The source code as a string
      virtual std::string to_string(const ASTResolver* resolver = nullptr) const = 0;

      /// Resolve the dependencies of a tree using an ASTResolver,
      /// which keeps track of parameters, dependent sub-trees, etc.
      /// \param resolver resolver object to use 
      virtual void resolve( ASTResolver& resolver ) = 0;
      /// virtual descructor
      virtual ~IExpression() = default;
      /// Evaluate the expression using the tree, 
      /// will generally be very slow but ocassionally useful for debugging. 
      virtual std::complex<double> operator()() const = 0;
      virtual Expression clone() const = 0; 
  };

  /// \ingroup ExpressionEngine cclass Expression 
  /// Wrapper for shared_ptrs to virtual expressions for use with overloaded operators, 
  /// i.e. what expression trees should be constructed from. 
  class Expression { 
    public:
    Expression();
    Expression( const double& value );
    Expression( const std::complex<double>& value );
    Expression( const std::shared_ptr<IExpression>& expression ) ;
    ~Expression() = default;
    std::string to_string(const ASTResolver* resolver = nullptr) const;
    IExpression* get() const;
    void resolve( ASTResolver& resolver );
    Expression operator+=( const Expression& other ) const;
    Expression operator*=( const Expression& other ) const;
    Expression operator-=( const Expression& other ) const;
    Expression operator/=( const Expression& other ) const;
    Expression operator-() const;
    Expression clone() const { return m_expression->clone(); } 
    std::complex<double> operator()() const; 

  private:
    std::shared_ptr<IExpression> m_expression;
  };

  /// \ingroup ExpressionEngine class Constant
  /// Class to contain a complex valued complex 
  class Constant : public IExpression {
    public:
    Constant( const int& real, const int& imag=0 ) : m_value(real,imag) {}
    Constant( const double& real, const double& imag=0) : m_value( real, imag ) {}
    Constant( const std::complex<double>& cmplx ) : m_value( cmplx) {}
    std::string to_string(const ASTResolver* resolver = nullptr) const override;
    void resolve( ASTResolver& resolver ) override;
    operator Expression() const;
    std::complex<double> operator()() const override { return m_value; }
    Expression clone() const override { return Constant(m_value); }
    private:
    std::complex<double> m_value;
  };

  /// \ingroup ExpressionEngine class Parameter
  /// \brief Free parameter for expression 
  ///
  /// Class to contain a function parameter. In general, this will include global function parameters 
  /// such as masses and widths, as well as the event buffer (i.e. kinematic quantities), which should be stored in 
  /// two separated parameter packs. There is also limited support for handling more complex function parameters to functions, 
  /// such as cache states, but this currently requires manually specifying the argument ordering. 
  struct Parameter : public IExpression {
    Parameter( const std::string& name="", const double& defaultValue = 0, const bool& resolved = false,
        const unsigned int& fromArg = 0 );
    std::string to_string(const ASTResolver* resolver = nullptr) const override;
    void resolve( ASTResolver& resolver ) override;
    operator Expression() const ;
    std::complex<double> operator()() const override { return std::complex<double>( m_defaultValue, 0 ); }
    std::string name() const { return m_name; }
    Expression clone() const override;
    std::string  m_name;
    bool         m_resolved;
    bool         m_compileTimeConstant;
    unsigned int m_fromArg;
    unsigned int m_address;
    double       m_defaultValue;
  };

  /// \ingroup ExpressionEngine class Ternary 
  /// Evaluates the ternary operator, i.e. 
  /// \code{.cpp}
  /// return a ? b : c 
  /// \endcode 
  struct Ternary : public IExpression {
    Ternary( const Expression& cond, const Expression& v1, const Expression& v2 );
    std::string to_string(const ASTResolver* resolver=nullptr) const override;
    void resolve( ASTResolver& resolver ) override;
    operator Expression() const ;
    std::complex<double> operator()() const override { return std::real(m_cond()) ? m_v1() : m_v2(); }
    Expression clone() const override { return Ternary( m_cond.clone(), m_v1.clone(), m_v2.clone() ) ; } 
    Expression m_cond;
    Expression m_v1;
    Expression m_v2;
  };

  /// \ingroup ExpressionEngine class SubTree
  struct SubTree : public IExpression { 
    SubTree( const Expression& other ) ;
    std::string to_string(const ASTResolver* resolver=nullptr) const override ;
    void resolve( ASTResolver& resolver ) override;
    operator Expression() const ;
    std::complex<double> operator()() const override { return m_expression(); }
    Expression clone() const override { return SubTree( m_expression.clone() ); }
    
    uint64_t key() const; 
    void setKey( const uint64_t& k ) ; 
    Expression  m_expression;
    std::string m_name;
    uint64_t    m_key; 
  };

  struct Function : public IExpression {
    Function( const std::string& name, const std::vector<Expression>& args ) ;
    std::string to_string(const ASTResolver* resolver=nullptr) const override ;
    void resolve( ASTResolver& resolver ) override;
    operator Expression() const ;
    std::complex<double> operator()() const override { return 0; }  
    Expression clone() const override { 
      std::vector<Expression> cloned_args; 
      for( auto& arg : m_args ) cloned_args.push_back( arg.clone() );
      return Function( m_name, cloned_args );
    } 
    std::string m_name;
    std::vector<Expression> m_args;
  };

  /// \ingroup ExpressionEngine class IBinaryExpression
  ///  Base class for binary expressions, i.e. those that take a pair of arguments (such as \f$+,-,\times,/\f$)
  class IBinaryExpression : public IExpression {
    public:
    IBinaryExpression( const Expression& l, const Expression& r ) : lval( l ), rval( r ){};
    void resolve( ASTResolver& resolver ) override;
    std::complex<double> operator()() const override = 0;
    Expression l() const { return lval ; }
    Expression r() const { return rval ; }
    protected:
    Expression lval;
    Expression rval;
  };

  /// \ingroup ExpressionEngine class Sum 
  /// \brief Binary expression that returns \f$l+r\f$
  DECLARE_BINARY_OPERATOR( Sum );

  /// \ingroup ExpressionEngine class Sub
  /// \brief Binary expression that returns \f$l-r\f$
  DECLARE_BINARY_OPERATOR( Sub );
  
  /// \ingroup ExpressionEngine class Sub
  /// \brief Binary expression that returns \f$l\times r\f$
  DECLARE_BINARY_OPERATOR( Product );
  
  /// \ingroup ExpressionEngine class Divide
  /// \brief Binary expression that returns \f$l / r\f$
  DECLARE_BINARY_OPERATOR( Divide );
  
  /// \ingroup ExpressionEngine class Pow
  /// \brief Binary expression that returns \f$ l^{r} \f$
  DECLARE_BINARY_OPERATOR( Pow );
  
  /// \ingroup ExpressionEngine class Fmod 
  /// \brief Binary expression that returns the fractional part of \f$ l \mathrm{mod} r \f$
  DECLARE_BINARY_OPERATOR( Fmod );
  
  /// \ingroup ExpressionEngine class LessThan
  /// \brief Binary expression that returns \f$l < r\f$
  DECLARE_BINARY_OPERATOR( LessThan );
  
  /// \ingroup ExpressionEngine class GreaterThan
  /// \brief Binary expression that returns \f$l > r\f$
  DECLARE_BINARY_OPERATOR( GreaterThan );
  
  /// \ingroup ExpressionEngine class And
  /// \brief Binary expression that returns \f$l \wedge r\f$
  DECLARE_BINARY_OPERATOR( And );
  DECLARE_BINARY_OPERATOR( Equal );

  /// \ingroup ExpressionEngine class IUnaryExpression
  ///  Base class for unary expressions, i.e. those that take a single argument.
  class IUnaryExpression : public IExpression {
    public:
    IUnaryExpression( const Expression& other ) : m_expression( other ){};
    void resolve( ASTResolver& resolver ) override;
    std::complex<double> operator()() const override = 0;
    virtual Expression d() const = 0; 
    Expression arg() const { return m_expression ;} 
    protected:
    Expression m_expression;
  };
  /// \ingroup ExpressionEngine struct Log
  /// \brief Unary expression that returns \f$\log(x)\f$
  DECLARE_UNARY_OPERATOR( Log );

  /// \ingroup ExpressionEngine struct Exp
  /// \brief Unary expression that returns \f$e^x\f$
  DECLARE_UNARY_OPERATOR( Exp );

  /// \ingroup ExpressionEngine struct Sqrt
  /// \brief Unary expression that returns \f$\sqrt{x}\f$
  DECLARE_UNARY_OPERATOR( Sqrt );

  /// \ingroup ExpressionEngine struct Abs 
  /// \brief Unary expression that returns \f$|z|\f$
  DECLARE_UNARY_OPERATOR( Abs );

  /// \ingroup ExpressionEngine struct Norm
  /// \brief Unary expression that returns \f$|z|^2\f$
  DECLARE_UNARY_OPERATOR( Norm );

  /// \ingroup ExpressionEngine struct Conj
  /// \brief Unary expression that returns \f$z^{*}\f$
  DECLARE_UNARY_OPERATOR( Conj );

  /// \ingroup ExpressionEngine struct Real
  /// \brief Unary expression that returns the real part of \f$z\f$
  DECLARE_UNARY_OPERATOR( Real );

  /// \ingroup ExpressionEngine struct Imag
  /// \brief Unary expression that returns the imaginary part of \f$z\f$
  DECLARE_UNARY_OPERATOR( Imag );

  /// \ingroup ExpressionEngine struct Sin
  /// \brief Unary expression that returns \f$\sin(z)\f$
  DECLARE_UNARY_OPERATOR( Sin );

  /// \ingroup ExpressionEngine struct Cos
  /// \brief Unary expression that returns \f$\cos(z)\f$
  DECLARE_UNARY_OPERATOR( Cos );

  /// \ingroup ExpressionEngine struct Tan
  /// \brief Unary expression that returns \f$\tan(z)\f$
  DECLARE_UNARY_OPERATOR( Tan );

  /// \ingroup ExpressionEngine struct ASin
  /// \brief Unary expression that returns \f$\sin^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ASin );

  /// \ingroup ExpressionEngine struct ACos
  /// \brief Unary expression that returns \f$\cos^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ACos );

  /// \ingroup ExpressionEngine struct ATan
  /// \brief Unary expression that returns \f$\tan^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ATan );

  Expression operator<( const Expression& A, const Expression& B );
  Expression operator>( const Expression& A, const Expression& B );

  Expression operator+( const Expression& A, const Expression& B );
  Expression operator-( const Expression& A, const Expression& B );
  Expression operator*( const Expression& A, const Expression& B );
  Expression operator/( const Expression& A, const Expression& B );

  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator+( const Expression& A, const T& B ){ return A + Constant(B); }
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator-( const Expression& A, const T& B ){ return A - Constant(B); }
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator*( const Expression& A, const T& B ){ return A * Constant(B); } 
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator/( const Expression& A, const T& B ){ return A / Constant(B); }

  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator+( const T& A, const Expression& B ){ return Constant(A) + B; }
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator-( const T& A, const Expression& B ){ return Constant(A) - B; }
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator*( const T& A, const Expression& B ){ return Constant(A) * B; }
  template < class T, typename std::enable_if_t< hasConstructor<Constant, T>() > >
  Expression operator/( const T& A, const Expression& B ){ return Constant(A) / B; }

  Expression operator&&( const Expression& A, const Expression& B );
  Expression operator==( const Expression& A, const Expression& B );
  Expression operator==( const Expression& A, const double& B );
  Expression operator==( const double& A, const Expression& B );

  std::ostream& operator<<( std::ostream& os, const Expression& expression ) ;
  namespace fcn {  
    Expression sqrt( const Expression& expression ); 
    Expression safe_sqrt( const Expression& expression );
    Expression complex_sqrt( const Expression& expression );
    Expression cos( const Expression& expression );
    Expression sin( const Expression& expression );
    Expression abs( const Expression& expression );
    Expression pow( const Expression& expression, const Expression& co );
    Expression norm( const Expression& expression );
    Expression conj( const Expression& expression );
    Expression exp( const Expression& expression );
    Expression log( const Expression& expression );
  }

  template < class T > bool is( const Expression& expression ){
    return dynamic_cast< const T* >( expression.get() ) != nullptr;
  }
  template < class T > T cast( const Expression& expression ){
    return *dynamic_cast< const T*>(expression.get() );
  }
  
  Expression make_cse( const Expression& A , bool simplify = false);

} // namespace AmpGen

#endif
