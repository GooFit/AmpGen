#include "AmpGen/CacheTransfer.h"

#include <iostream>
#include <string>

#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"

using namespace AmpGen;

CacheTransfer::CacheTransfer() = default; 

CacheTransfer::CacheTransfer( const size_t& address, const std::string& name, const double& value, const size_t& size ) : 
  m_address(address), 
  m_size(size),
  m_value(value),
  m_name(name)
{
} 

void CacheTransfer::transfer( CompiledExpressionBase* destination )
{ 
  destination->setExternal(m_value, m_address);
}

void CacheTransfer::print() const 
{ 
  INFO( m_address << " " << m_value << " " << m_name ) ; 
}

void ParameterTransfer::transfer( CompiledExpressionBase* destination )
{
  destination->setExternal( m_source->mean(), m_address );
}

ParameterTransfer::ParameterTransfer(const size_t& address, const std::string& name, MinuitParameter* source )
  : CacheTransfer(address, name, source->mean(), 1), 
    m_source( source )
{
}

void ParameterTransfer::print() const 
{ 
  INFO( "Source: " << m_source->name() << " address = " << m_address << " value = " << m_source->mean() ); 
}
