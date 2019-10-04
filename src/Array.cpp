#include "AmpGen/Array.h"

#include <complex>
#include <vector>

#include "AmpGen/ASTResolver.h"
#include "AmpGen/CacheTransfer.h"

using namespace AmpGen; 

Array::Array( const Expression& top, 
    const size_t&     size,
    const Expression& address)
  : m_top(top),
  m_address(address),
  m_size(size)
{
}

std::string Array::to_string(const ASTResolver* resolver) const {
  auto head   = m_top.to_string(resolver);
  if( is<Constant>(m_address) ) 
    return head+"["+ std::to_string(int(std::real(m_address()))) +"]";
  auto offset = m_address.to_string(resolver);
  auto pos = head.find_last_of("]");
  if( pos != std::string::npos ){
    auto st1 = head.substr(0,pos);
    auto st2 = head.substr(pos,head.size() - pos );
    return st1 + "+int("+offset+")" + st2;
  }
  else return head +"[ int("+offset+")]";
}

void Array::resolve( ASTResolver& resolver ) const
{
  m_top.resolve( resolver );
  m_address.resolve( resolver );
}

Array::operator Expression() { return Expression( std::make_shared<Array>( *this ) ); }
complex_t Array::operator()() const { return 0; }

Expression Array::operator[]( const Expression& address ) const { 
  return Array( m_top, m_size, address );
}
