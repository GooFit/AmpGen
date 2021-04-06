#ifndef AMPGEN_CACHETRANSFER_H
#define AMPGEN_CACHETRANSFER_H

#include <vector>
#include <cstddef>
#include <string>
#include <functional>

namespace AmpGen
{
  class CompiledExpressionBase;
  class MinuitParameter;
  class LambdaExpression; 

  class CacheTransfer 
  {
    public:
      CacheTransfer();
      CacheTransfer( const size_t& address, const std::string& name, const double& value=0, const size_t& size=1);
      virtual ~CacheTransfer() = default;

      size_t address() const { return m_address ; }

      virtual void transfer( CompiledExpressionBase* destination ); 
      virtual void print()     const;
      virtual size_t size()    const { return m_size ; }  
      virtual std::string name() const { return m_name; } 
    protected: 
      size_t       m_address = {0};
      size_t       m_size    = {0}; 
      double       m_value   = {0};
      std::string  m_name    = {""}; 
  };

  class ParameterTransfer : public CacheTransfer
  {
    public:
      ParameterTransfer( const size_t& address, const std::string& name, MinuitParameter* source );
      virtual ~ParameterTransfer() = default;

      size_t size()    const override { return 1 ; }  

      void transfer( CompiledExpressionBase* destination ) override;
      void print()     const override;
      virtual std::string name() const override; 

    protected:
      MinuitParameter* m_source = {nullptr};
  };
  class LambdaTransfer : public CacheTransfer 
  {
    public: 
      LambdaTransfer( const size_t& address, const std::string& name, const LambdaExpression* source );

      size_t size()    const override { return 1 ; }  

      void transfer( CompiledExpressionBase* destination ) override;
      void print()     const override;

      std::function<double(void)> m_function;   
  };

} // namespace AmpGen

#endif
