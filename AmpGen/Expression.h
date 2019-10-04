#ifndef AMPGEN_EXPRESSION_H
#define AMPGEN_EXPRESSION_H

/** @defgroup ExpressionEngine Expressions 
    @brief are a set of classes for constructing mathematical functions from primitives.
   
    AmpGen Expressions are a set of classes for constructing functions from primatives 
    such as the binary operators \f$(+,-,/,\times)\f$, and elementary functions such as trigonometric 
    functions and sqrts. 
    These are then compiled to conventional C++ at runtime, to allow fast and flexible 
    evalulation of amplitudes. 
    Expression trees can be built using overloaded operators, so 

    \code{.cpp}
    Expression x = a+b; 
    \endcode

    will given an expression that can calculate the sum of expressions a and b.
    There are constructors for expressions from constants and corresponding overloads of operators, so for example
    
    \code{.cpp}
    Expression x = a+5;
    \endcode

    generates a tree that adds '5' to the the tree 'a'.  
    Assorted optimisations will generally be invoked on the AST if constants are included, such as reduction of products and sums of constants.
    
    There is limited functionality for expressions to be evaluated in-situ to the default values associated with their parameters, which can useful for debugging purposes.
    The more typical usage is to use an expression to generate a CompiledExpression that compiles the expression and wraps the corresponding shared library with a functor. 

    For unary functions, it is recommend to use the functions contained in the fcn:: namespace that have the same names as those in std::math. For example, 
    \code{.cpp}
    Expression y = fcn::sin(x);
    \endcode

    Use of the unary classes directly can cause unexpected behaviour if used in conjunction with the auto keyword, and does not include some optimisations as readily such as constant folding. 
    */

/// @ingroup ExpressionEngine macro ADD_DEBUG
/// Make a (named) debugging expression and add to a set of DebugSymbols.
#define ADD_DEBUG( X, Y )                                               \
  if ( Y != 0 ) Y->push_back( DebugSymbol( std::string( #X ), X ) );    \

/// @ingroup ExpressionEngine macro DEFINE_CAST 
/// Define a cast from an expression implementation to the shared_ptr wrapper. 
#define DEFINE_CAST( X )                                    \
  X::operator Expression() const                            \
{ return Expression( std::make_shared<X>( *this ) ); }  

#define DEFINE_BINARY_OPERATOR( X ) \
  X::X( const AmpGen::Expression& l, const AmpGen::Expression& r ) : IBinaryExpression(l,r ) {} \
  X::operator Expression() const { return Expression( std::make_shared<X>(*this) ) ; } 

#define DEFINE_UNARY_OPERATOR( X, F ) \
  X::X( const AmpGen::Expression& expression) : IUnaryExpression(expression) {} \
  X::operator Expression() const { return Expression( std::make_shared<X>(*this) ) ; } \
  complex_t X::operator()() const { return F( m_expression() ); } \
  std::string X::to_string(const ASTResolver* resolver) const { return std::string(#F)+"("+ m_expression.to_string(resolver)+")";}

/// @ingroup ExpressionEngine macro DECLARE_UNARY_OPERATOR
/// Macro to declare a unary operator, \ref ExpressionEngine "see IUnaryExpression"
#define DECLARE_UNARY_OPERATOR( X )                         \
  class X : public IUnaryExpression {                       \
    public:                                                 \
    explicit X( const Expression& other );                  \
    virtual std::string to_string(const ASTResolver* resolver=nullptr) const override;         \
    virtual Expression d() const override;                  \
    operator Expression() const;                            \
    virtual complex_t operator()() const override;          \
  }

/// @ingroup ExpressionEngine macro DECLARE_BINARY_OPERATOR
/// Macro to declare a binary operator, \ref ExpressionEngine "see IBinaryExpression"
#define DECLARE_BINARY_OPERATOR( X )                        \
  class X : public IBinaryExpression {                      \
    public:                                                 \
    X( const Expression& l, const Expression& r );          \
    virtual std::string to_string(const ASTResolver* resolver=nullptr) const override ;        \
    operator Expression() const ;                           \
    virtual complex_t operator()() const override;          \
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
#include "AmpGen/Types.h"

namespace AmpGen
{
  class ASTResolver;
  class Expression;
  class Parameter;

  typedef std::pair<std::string, Expression> DebugSymbol;
  typedef std::vector<DebugSymbol> DebugSymbols;

  /** @ingroup ExpressionEngine class IExpression 
      @brief Virtual base class for other expression tree components. <br>
      Implementations must permit the following operations on nodes. */
  class IExpression {
    public:
      /// Called to convert the Expression tree into source code.
      /// \return The source code as a string
      virtual std::string to_string(const ASTResolver* resolver=nullptr) const = 0;
      /// Resolve the dependencies of a tree using an ASTResolver,
      /// which keeps track of parameters, dependent sub-trees, etc.
      /// \param resolver resolver object to use 
      virtual void resolve( ASTResolver& resolver ) const = 0;
      /// virtual descructor
      virtual ~IExpression() = default;
      /// Evaluate the expression using the tree, 
      /// will generally be very slow but ocassionally useful for debugging. 
      virtual complex_t operator()() const = 0;
  };

  /** @ingroup ExpressionEngine class Expression 
      @brief Wrapper class for shared_ptrs to virtual expressions for use in conjunction with operators 
      to build expression trees. */
  class Expression { 
    public:
    Expression();
    Expression( const real_t& value );
    Expression( const complex_t& value );
    Expression( const std::shared_ptr<IExpression>& expression ) ;
    ~Expression() = default;
    std::string to_string(const ASTResolver* resolver=nullptr) const;
    IExpression* get() const;
    void resolve( ASTResolver& resolver ) const;
    Expression operator+=( const Expression& other ) const;
    Expression operator*=( const Expression& other ) const;
    Expression operator-=( const Expression& other ) const;
    Expression operator/=( const Expression& other ) const;
    Expression operator-() const;
    complex_t operator()() const; 

  private:
    std::shared_ptr<IExpression> m_expression;
  };

  /** @ingroup ExpressionEngine class Constant
      @brief Class to contain a constant (which can contain a complex value) */
   class Constant : public IExpression {
    public:
    template <typename T1, 
              typename T2,
              typename = std::enable_if_t< std::is_constructible< complex_t, T1, T2 >::value > >
    Constant( const T1& real, const T2& imag=0 ) : m_value(real,imag) {}

    Constant( const complex_t& value ) : m_value(value) {}
    std::string to_string(const ASTResolver* resolver=nullptr) const override;
    void resolve( ASTResolver& resolver ) const override;
    operator Expression() const;
    complex_t operator()() const override { return m_value; }
    private:
    complex_t m_value;
  };

  /** @ingroup ExpressionEngine class Parameter
      @brief Free parameter for expression 
  
   Class to contain a parameter. In general, this will include global function parameters 
   such as masses and widths, as well as the event buffer (i.e. kinematic quantities), which should be stored in 
   two separated parameter packs. There is also limited support for handling more complex function parameters to functions, 
   such as cache states, but this currently requires manually specifying the argument ordering. */
  class Parameter : public IExpression {
    public:
    Parameter(const std::string& name     = "", 
              const double&  defaultValue = 0 , 
              const bool&        resolved = false);
    std::string to_string(const ASTResolver* resolver = nullptr) const override;
    void resolve( ASTResolver& resolver ) const override;
    operator Expression() const;
    complex_t operator()() const override { return complex_t( m_defaultValue, 0 ); }
    std::string name() const { return m_name; }
    const double& defaultValue() const { return m_defaultValue ; }
    double& defaultValue() { return m_defaultValue ; }
    bool isResolved() const { return m_resolved ;}
    private:
    std::string  m_name;
    double       m_defaultValue;
    bool         m_resolved; 
  };

  /** @ingroup ExpressionEngine class Ternary 
      @brief Evaluates the ternary operator.

      The equivalent c++ code would be
      \code{.cpp}
      return a ? b : c 
      \endcode */
  struct Ternary : public IExpression {
    Ternary( const Expression& cond, const Expression& v1, const Expression& v2 );
    std::string to_string(const ASTResolver* resolver = nullptr ) const override;
    void resolve( ASTResolver& resolver ) const override;
    operator Expression() const ;
    complex_t operator()() const override { return std::real(m_cond()) ? m_v1() : m_v2(); }
    Expression m_cond;
    Expression m_v1;
    Expression m_v2;
  };

  /// @ingroup ExpressionEngine class SubTree
  struct SubTree : public IExpression { 
    SubTree( const Expression& other ) ;
    std::string to_string(const ASTResolver* resolver = nullptr ) const override ;
    void resolve( ASTResolver& resolver ) const override;
    operator Expression() const ;
    complex_t operator()() const override { return m_expression(); }
    uint64_t key() const;
    void setKey( const size_t& new_key ); 
    Expression  m_expression;
    uint64_t    m_key; 
  };

  struct Function : public IExpression {
    Function( const std::string& name, const std::vector<Expression>& args ) ;
    std::string to_string(const ASTResolver* resolver = nullptr ) const override ;
    void resolve( ASTResolver& resolver ) const override;
    operator Expression() const ;
    complex_t operator()() const override { return 0; }  
    std::string m_name;
    std::vector<Expression> m_args;
  };

  /// @ingroup ExpressionEngine class IBinaryExpression
  ///  Base class for binary expressions, i.e. those that take a pair of arguments (such as \f$+,-,\times,/\f$)
  class IBinaryExpression : public IExpression {
    public:
    IBinaryExpression( const Expression& l, const Expression& r ) : lval( l ), rval( r ){};
    void resolve( ASTResolver& resolver ) const override;
    complex_t operator()() const override = 0;
    Expression l() const { return lval ; }
    Expression r() const { return rval ; }
    protected:
    Expression lval;
    Expression rval;
  };

  /// @ingroup ExpressionEngine class Sum 
  /// @brief Binary expression that returns \f$l+r\f$
  DECLARE_BINARY_OPERATOR( Sum );

  /// @ingroup ExpressionEngine class Sub
  /// @brief Binary expression that returns \f$l-r\f$
  DECLARE_BINARY_OPERATOR( Sub );
  
  /// @ingroup ExpressionEngine class Sub
  /// @brief Binary expression that returns \f$l\times r\f$
  DECLARE_BINARY_OPERATOR( Product );
  
  /// @ingroup ExpressionEngine class Divide
  /// @brief Binary expression that returns \f$l / r\f$
  DECLARE_BINARY_OPERATOR( Divide );
  
  /// @ingroup ExpressionEngine class Pow
  /// @brief Binary expression that returns \f$ l^{r} \f$
  DECLARE_BINARY_OPERATOR( Pow );
  
  /// @ingroup ExpressionEngine class Fmod 
  /// @brief Binary expression that returns the fractional part of \f$ l \mathrm{mod} r \f$
  DECLARE_BINARY_OPERATOR( Fmod );
  
  /// @ingroup ExpressionEngine class LessThan
  /// @brief Binary expression that returns \f$l < r\f$
  DECLARE_BINARY_OPERATOR( LessThan );
  
  /// @ingroup ExpressionEngine class GreaterThan
  /// @brief Binary expression that returns \f$l > r\f$
  DECLARE_BINARY_OPERATOR( GreaterThan );
  
  /// @ingroup ExpressionEngine class And
  /// @brief Binary expression that returns \f$l \wedge r\f$
  DECLARE_BINARY_OPERATOR( And );
  DECLARE_BINARY_OPERATOR( Equal );

  DECLARE_BINARY_OPERATOR( ATan2 );
  /// @ingroup ExpressionEngine class IUnaryExpression
  ///  Base class for unary expressions, i.e. those that take a single argument.
  class IUnaryExpression : public IExpression {
    public:
    IUnaryExpression( const Expression& other ) : m_expression( other ){};
    void resolve( ASTResolver& resolver ) const override;
    complex_t operator()() const override = 0;
    virtual Expression d() const = 0; 
    Expression arg() const { return m_expression ;} 
    protected:
    Expression m_expression;
  };
  /// @ingroup ExpressionEngine struct Log
  /// @brief Unary expression that returns \f$\log(x)\f$
  DECLARE_UNARY_OPERATOR( Log );

  /// @ingroup ExpressionEngine struct Exp
  /// @brief Unary expression that returns \f$e^x\f$
  DECLARE_UNARY_OPERATOR( Exp );

  /// @ingroup ExpressionEngine struct Sqrt
  /// @brief Unary expression that returns \f$\sqrt{x}\f$
  DECLARE_UNARY_OPERATOR( Sqrt );

  /// @ingroup ExpressionEngine struct Abs 
  /// @brief Unary expression that returns \f$|z|\f$
  DECLARE_UNARY_OPERATOR( Abs );

  /// @ingroup ExpressionEngine struct Norm
  /// @brief Unary expression that returns \f$|z|^2\f$
  DECLARE_UNARY_OPERATOR( Norm );

  /// @ingroup ExpressionEngine struct Conj
  /// @brief Unary expression that returns \f$z^{*}\f$
  DECLARE_UNARY_OPERATOR( Conj );

  /// @ingroup ExpressionEngine struct Real
  /// @brief Unary expression that returns the real part of \f$z\f$
  DECLARE_UNARY_OPERATOR( Real );

  /// @ingroup ExpressionEngine struct Imag
  /// @brief Unary expression that returns the imaginary part of \f$z\f$
  DECLARE_UNARY_OPERATOR( Imag );

  /// @ingroup ExpressionEngine struct Sin
  /// @brief Unary expression that returns \f$\sin(z)\f$
  DECLARE_UNARY_OPERATOR( Sin );

  /// @ingroup ExpressionEngine struct Cos
  /// @brief Unary expression that returns \f$\cos(z)\f$
  DECLARE_UNARY_OPERATOR( Cos );

  /// @ingroup ExpressionEngine struct Tan
  /// @brief Unary expression that returns \f$\tan(z)\f$
  DECLARE_UNARY_OPERATOR( Tan );

  /// @ingroup ExpressionEngine struct ASin
  /// @brief Unary expression that returns \f$\sin^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ASin );

  /// @ingroup ExpressionEngine struct ACos
  /// @brief Unary expression that returns \f$\cos^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ACos );
  /// @ingroup ExpressionEngine struct ATan
  /// @brief Unary expression that returns \f$\tan^{-1}(z)\f$
  DECLARE_UNARY_OPERATOR( ATan );
  
  /// @ingroup ExpressionEngine struct ISqrt
  /// @brief Unary expression that uses inverse sqrt (faster than dividing by sqrt) 
  DECLARE_UNARY_OPERATOR( ISqrt );
  
  /// @ingroup ExpressionEngine struct LGamma
  /// @brief Unary expression that returns \f$\log\left|\Gamma(x)\right|\f$
  DECLARE_UNARY_OPERATOR( LGamma );

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
    Expression isqrt( const Expression& expression );
    Expression cos( const Expression& expression );
    Expression sin( const Expression& expression );
    Expression tan( const Expression& expression );
    Expression abs( const Expression& expression );
    Expression acos( const Expression& expression );
    Expression asin( const Expression& expression );
    Expression atan( const Expression& expression );
    Expression pow( const Expression& expression, const Expression& co );
    Expression fpow( const Expression& expression, const int& n );
    Expression norm( const Expression& expression );
    Expression conj( const Expression& expression );
    Expression exp( const Expression& expression );
    Expression log( const Expression& expression );
    Expression atan2( const Expression& y, const Expression& x);
  }

  template < class T > bool is( const Expression& expression ){
    return dynamic_cast< const T* >( expression.get() ) != nullptr;
  }
  template < class T > T cast( const Expression& expression ){
    return *static_cast< const T*>(expression.get() );
  }
  
  Expression make_cse( const Expression& A , bool simplify = false);

} // namespace AmpGen

#endif
