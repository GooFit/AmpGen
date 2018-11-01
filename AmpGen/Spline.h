#ifndef AMPGEN_SPLINE_H
#define AMPGEN_SPLINE_H

#include "AmpGen/Array.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CacheTransfer.h"

#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include "TMatrixD.h"

namespace AmpGen{
  class SplineTransfer : public CacheTransfer
  {
  private:
    TMatrixD m_transferMatrix;
    std::vector<AmpGen::MinuitParameter*> m_parameters;
    unsigned int m_nKnots;
    double       m_min;
    double       m_max;
    unsigned int m_address;

  public:
    SplineTransfer();
    SplineTransfer( const SplineTransfer& other );
    SplineTransfer( const unsigned int& address, const unsigned int& N, const double& min, const double& max );
    void transfer( CompiledExpressionBase* destination ) override;

    bool isConfigured();

    void set( const unsigned int& N, AmpGen::MinuitParameter* f );
    void set( const unsigned int& N, const double& value );

    void setAddress( const unsigned int& address );
    void print() const override;
    virtual unsigned int address() const { return m_address ; } 
    virtual unsigned int size() const { return 2*m_nKnots ; }  
  };
  struct Spline {
    Array       m_points; 
    std::string m_name;
    size_t      m_nKnots;
    double      m_min;
    double      m_max;
    std::vector<real_t> m_values; 
    Spline(const std::string& name, const size_t& nBins, const double& min, const double& max, const std::vector<double>& values = {} ) ;
    Expression eval(const Expression& x , DebugSymbols* db=nullptr);
    void resolve( ASTResolver& resolver );
    Expression operator()( const Expression& x, DebugSymbols* db=nullptr );
    void set( const std::vector<real_t>& values );
    Spline clone() const { return Spline(m_name,m_nKnots,m_min,m_max) ; } 
  };
  
  struct SplineExpression : public IExpression { 
    std::shared_ptr<Spline> m_parent; 
    Expression                   m_x; 
    Expression                   m_eval;
    SplineExpression( const Spline& parent, const Expression& x , DebugSymbols* db = nullptr) ; 
    void resolve( ASTResolver& resolver ) override ;
    std::string to_string(const ASTResolver* resolver=nullptr) const { return m_eval.to_string(resolver) ; } 
    operator Expression() { return Expression( std::make_shared<SplineExpression>( *this ) ); }
    complex_t operator()() const override { return 0; }
    Expression clone() const { return SplineExpression( m_parent->clone(), m_x.clone() ) ; }

  };
  Expression getSpline( const std::string& name, const AmpGen::Expression& x, const std::string& arrayName,
                              AmpGen::DebugSymbols* dbexpressions = nullptr, const bool& continueSpline = false );

}

#endif 
