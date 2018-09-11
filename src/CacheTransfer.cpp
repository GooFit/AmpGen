#include "AmpGen/CacheTransfer.h"

#include <TMatrixTUtils.h>
#include <TVectorDfwd.h>
#include <TVectorT.h>
#include <ostream>

#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

using namespace AmpGen;

CacheTransfer::CacheTransfer( const unsigned int& address, const double& value, const size_t& size ) : 
  m_address(address), m_value(value), m_size(size) {} 

void CacheTransfer::transfer( CompiledExpressionBase* destination )
{ 
  destination->setExternal(m_value, m_address);
}
void CacheTransfer::print() const { INFO( m_address << " " << m_value ) ; }

void ParameterTransfer::transfer( CompiledExpressionBase* destination )
{
  destination->setExternal( m_source->mean(), m_address );
}

CacheTransfer::CacheTransfer() : m_address( 0 ), m_value( 0 ) {}

ParameterTransfer::ParameterTransfer( const unsigned int& address, MinuitParameter* source )
  : m_address( address ), m_source( source )
{
}

void ParameterTransfer::print() const { std::cout << this << " " << m_source->name() << " " << m_address << std::endl; INFO( "Source: " << m_source->name() << " address = " << m_address ); }
