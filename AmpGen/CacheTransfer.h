#ifndef AMPGEN_CACHETRANSFER_H
#define AMPGEN_CACHETRANSFER_H

#include <vector>
#include <cstddef>

namespace AmpGen
{
  class CompiledExpressionBase;
  class MinuitParameter;
  class CacheTransfer 
  {
    protected: 
      unsigned int m_address;
      double       m_value; 
      size_t       m_size; 
    public:
      CacheTransfer();
      CacheTransfer( const unsigned int& address, const double& value, const size_t& size=1);
      virtual ~CacheTransfer() = default;
      virtual void transfer( CompiledExpressionBase* destination ) ; 
      virtual void print() const ;
      virtual unsigned int address() const { return m_address ; }
      virtual unsigned int size() const { return m_size ; }  
  };

  class ParameterTransfer : public CacheTransfer
  {
  protected:
    unsigned int m_address;
    AmpGen::MinuitParameter* m_source;

  public:
    ParameterTransfer( const unsigned int& address, AmpGen::MinuitParameter* source );
    virtual ~ParameterTransfer() = default;
    virtual void transfer( CompiledExpressionBase* destination );
    virtual void print() const;
    virtual unsigned int address() const { return m_address ; } 
    virtual unsigned int size() const { return 1 ; }  
  };

} // namespace AmpGen

#endif
