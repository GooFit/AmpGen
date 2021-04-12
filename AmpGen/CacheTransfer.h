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
    public:
      CacheTransfer();
      CacheTransfer( const size_t& address, const double& value=0, const size_t& size=1);
      virtual ~CacheTransfer() = default;
      
      size_t address() const { return m_address ; }
      
      virtual void transfer( CompiledExpressionBase* destination ); 
      virtual void print()     const;
      virtual size_t size()    const { return m_size ; }  
    
    protected: 
      size_t       m_address = {0};
      size_t       m_size    = {0}; 
      double       m_value   = {0}; 
  };

  class ParameterTransfer : public CacheTransfer
  {
  public:
    ParameterTransfer( const size_t& address, AmpGen::MinuitParameter* source );
    virtual ~ParameterTransfer() = default;
    
    size_t size()    const override { return 1 ; }  
    
    void transfer( CompiledExpressionBase* destination ) override;
    void print()     const override;
  
  protected:
    AmpGen::MinuitParameter* m_source = {nullptr};
  };

} // namespace AmpGen

#endif
