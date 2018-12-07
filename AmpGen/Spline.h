#ifndef AMPGEN_SPLINE_H
#define AMPGEN_SPLINE_H

#include <memory.h>
#include <stddef.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Array.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CacheTransfer.h"
#include "AmpGen/Types.h"

#include "TMatrixD.h"

namespace AmpGen{
  class ASTResolver;
  class CompiledExpressionBase;
  class MinuitParameter;

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
  
  struct Spline : public IExpression { 
    Array                        m_points; 
    std::string                  m_name;
    size_t                       m_nKnots;
    double                       m_min;
    double                       m_max;
    Expression                   m_x; 
    Expression                   m_eval; 
    Spline( const std::string& name, 
            const size_t& nKnots, 
            const double& min, 
            const double& max );
 
    Spline( const Spline& spline, const Expression& x );
    void resolve( ASTResolver& resolver ) override ;
    std::string to_string(const ASTResolver* resolver=nullptr) const ;
    operator Expression() ;
    complex_t operator()() const override ;
    Expression clone() const ;
    Expression operator()( const Expression& x ); 
    Expression eval() const ;
  };
  Expression getSpline( const std::string& name, const AmpGen::Expression& x, const std::string& arrayName,
                              AmpGen::DebugSymbols* dbexpressions = nullptr, const bool& continueSpline = false );
}

#endif 
