#ifndef AMPGEN_ARRAY_H
#define AMPGEN_ARRAY_H

#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpressionBase.h"

#include <vector>
#include <cstddef>

namespace AmpGen
{
  struct Array {
    Expression   m_top;
    size_t       m_address; 
    size_t       m_size; 
    Array( const Expression& expression, const size_t& size ) : 
      m_top(expression), 
      m_address(9999),
      m_size(size) {}
    Expression operator[]( const Expression& address );
    virtual ~Array() = default;
    void resolve( ASTResolver& resolver ) ;
    Array clone() const { 
      return Array( m_top.clone(), m_size );
    }
  };

  struct ArrayExpression : public IExpression {
    std::shared_ptr<Array> m_parent;
    Expression m_addressOffset;
    Expression clone() const {
      return ArrayExpression( m_parent->clone(), m_addressOffset.clone() );
    }

    ArrayExpression( const Array& parent, const Expression& address )
        : m_parent( std::make_shared<Array>( parent ) ), m_addressOffset( address )
    {
    }
    ArrayExpression( const std::shared_ptr<Array>& parent, const Expression& address )
        : m_parent(parent), m_addressOffset( address )
    {
    }
    std::string to_string(const ASTResolver* resolver=nullptr) const override ;

    void resolve( ASTResolver& resolver ) override
    {
      m_addressOffset.resolve( resolver );
      m_parent->resolve( resolver );
    }
    operator Expression() { return Expression( std::make_shared<ArrayExpression>( *this ) ); }
    std::complex<double> operator()() const override { return 0; }
  };
} // namespace AmpGen

#endif
