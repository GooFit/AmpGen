#include "AmpGen/Array.h"

using namespace AmpGen; 

AmpGen::Expression AmpGen::Array::operator[]( const AmpGen::Expression& expression )
{
  return AmpGen::Expression( AmpGen::ArrayExpression( *this, expression ) );
}

void Array::resolve( ASTResolver& resolver ) 
{
  double value = 0.0; 
  if( is<Parameter>(m_top) ){
    auto param = cast<Parameter>(m_top);
    param.resolve( resolver );
    if( ! param.m_resolved ) param.m_address = resolver.addCacheFunction<CacheTransfer>( param.m_name, value , m_size );
    m_address = param.m_address; 
  }
  else {
    m_address = 0;
    m_top.resolve( resolver );
  }
}

std::string ArrayExpression::to_string() const
{
  auto& top = m_parent->m_top;

  if( is<Parameter>(top) ){
    auto param = cast<Parameter>( top );
    std::string top_name = "x"+std::to_string(param.m_fromArg); 
    if( is<Constant>( m_addressOffset ) ){     
      return top_name + "["+ std::to_string( int( m_parent->m_address + std::real( m_addressOffset() ) ) ) +"]"; 
    }
    else 
      return top_name + "["+ std::to_string( m_parent->m_address ) + "+int("+m_addressOffset.to_string() +")]"; 
  }
  else {
    if( is<Constant>(m_addressOffset ) ) return top.to_string() + "["+std::to_string( int(std::real(m_addressOffset())) ) + "]";
    else return top.to_string() + "[int("+m_addressOffset.to_string() +")]";
  }
}
