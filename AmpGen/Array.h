#ifndef AMPGEN_ARRAY_H
#define AMPGEN_ARRAY_H

#include <memory.h>
#include <vector>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>

#include "AmpGen/Expression.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class ASTResolver;

  /** @ingroup ExpressionEngine class Array
      @brief Expression for a fixed size array of values.
      Expression for an array, i.e. a set of values with an index.
      Can be used to return an expression from the array, which is resolved at (second) compile time. 
      */  
  class Array : public IExpression 
  {
    public:
      Array( const Expression& top, const size_t& size, const Expression& address = 0 );
      std::string to_string(const ASTResolver* resolver=nullptr) const override ;
      void resolve( ASTResolver& resolver ) const override;
      operator Expression(); 
      complex_t operator()() const override;
      Expression operator[]( const Expression& address ) const;
      Expression top() const { return m_top ; } 

    private:
      Expression m_top; 
      Expression m_address;
      size_t     m_size; 
  };
} // namespace AmpGen

#endif
