#ifndef AMPGEN_SPLINE_H
#define AMPGEN_SPLINE_H

#include <memory.h>
#include <stddef.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gsl/gsl_matrix.h>

#include "AmpGen/Array.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CacheTransfer.h"
#include "AmpGen/Types.h"

namespace AmpGen{
  class ASTResolver;
  class CompiledExpressionBase;
  class MinuitParameter;

  class SplineTransfer : public CacheTransfer
  {
    public:
//      SplineTransfer();
      SplineTransfer( const size_t& address, const std::string& name, const unsigned int& N, const double& min, const double& max );
      void transfer( CompiledExpressionBase* destination ) override;
      bool isConfigured();
      void set( const unsigned int& N, MinuitParameter* f );
      void set( const unsigned int& N, const double& value );
      void print()  const override;
      size_t size() const override { return 2*m_nKnots; }
      
      SplineTransfer( const SplineTransfer&) = delete;
      SplineTransfer& operator=(const SplineTransfer&) = delete; 
      SplineTransfer& operator=(SplineTransfer&&) = delete; 
      SplineTransfer( SplineTransfer&&) = delete; 
      ~SplineTransfer(); 
    private:
      gsl_matrix*  m_transferMatrix = {nullptr}; 
      std::vector<MinuitParameter*> m_parameters;
      size_t       m_nKnots;
      double       m_min;
      double       m_max;
  };

  class Spline : public IExpression { 
    public:
      Spline( const std::string& name, 
          const size_t& nKnots, 
          const double& min, 
          const double& max );

      Spline( const Spline& spline, const Expression& x, DebugSymbols* db =nullptr );
      void resolve( ASTResolver& resolver ) const override ;
      std::string to_string(const ASTResolver* resolver=nullptr) const override;
      operator Expression() ;
      complex_t operator()() const override ;
      Expression operator()( const Expression& x, DebugSymbols* db); 
      Expression eval(DebugSymbols* db=nullptr) const ;

      Array                        m_points; 
      std::string                  m_name;
      size_t                       m_nKnots;
      double                       m_min;
      double                       m_max;
      Expression                   m_x; 
      Expression                   m_eval; 
  };
  Expression getSpline( const std::string& name, const Expression& x, const std::string& arrayName,
      DebugSymbols* dbexpressions = nullptr, const bool& continueSpline = false );
}

#endif 
