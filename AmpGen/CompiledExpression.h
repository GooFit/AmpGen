#ifndef AMPGEN_COMPILEDEXPRESSION_H
#define AMPGEN_COMPILEDEXPRESSION_H

#include "AmpGen/CacheTransfer.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/DynamicFCN.h"
#include "AmpGen/Expression.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"
#include "AmpGen/simd/utils.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/ArgumentPack.h"
#include <cxxabi.h>
#include <dlfcn.h>
#include <vector>
#include <map>

namespace AmpGen
{
  /* @class CompiledExpression
     @tparam ret_type The type that is returned this compiled expression,
     usually this is a std::complex<double>,
     but in principal support also exists for computing coupled channel propagators
     (i.e. returning array types) */
  namespace detail {
    template <typename T> struct size_of {       static constexpr unsigned value = sizeof(T); };
    template <>           struct size_of<void> { static constexpr unsigned value = 0; } ;
  }
  DECLARE_ARGUMENT(disableBatch, bool);
  DECLARE_ARGUMENT(includeParameters, bool); 
  DECLARE_ARGUMENT(includePythonBindings, bool); 

  template <typename ret_type, typename... arg_types> class CompiledExpression; 
  template <typename ret_type, typename... arg_types>
    class CompiledExpression<ret_type(arg_types...)> : public CompiledExpressionBase
    {

      private:
        DynamicFCN<ret_type( arg_types... )>                                  m_fcn;
        DynamicFCN<void( const size_t&, const size_t&, const size_t&, ret_type*, arg_types... )> m_batchFcn;
        DynamicFCN<std::vector<std::pair<std::string, complex_v>>(arg_types...)>      m_fdb;
        std::vector<real_t>  m_externals             = {};
        bool                 m_hasExternalsChanged   = {false};
      public:
        typedef ret_type return_type;
        unsigned m_outputSize = {0};

        template <typename... namedArgs> CompiledExpression( const Expression& expression, const std::string& name, const namedArgs&... args ) : 
          CompiledExpressionBase(expression, name) 
        {
          set(expression,name, args...); 
        }

        template <typename... namedArgs> void set( const Expression& expression, const std::string& name, const namedArgs&... args )  
        {
          m_obj = expression; 
          m_name = name; 
          const MinuitParameterSet* mps = nullptr; 
          auto process_argument = [this, &mps]( const auto& arg ) mutable
          {
            DEBUG( type_string(arg) ); 
            if constexpr( std::is_convertible<decltype(arg), DebugSymbols>::value  ){
              this->m_db = arg;
            }
            else if constexpr( std::is_convertible<decltype(arg), std::map<std::string, unsigned>>::value ) this->m_evtMap = arg;
            else if constexpr( std::is_convertible<decltype(arg), const MinuitParameterSet*>::value or 
                               std::is_convertible<decltype(arg), const AmpGen::MinuitParameterSet*>::value or
                               std::is_convertible<decltype(arg), MinuitParameterSet*>::value ){
              mps = arg;
            }
            else if constexpr( std::is_convertible<decltype(arg), MinuitParameterSet>::value ) mps = &arg;
            
            else if constexpr( std::is_convertible<decltype(arg), disableBatch>::value ) {
              DEBUG("Disabling bulk evaluation: did you do this on purpose?");
              m_disableBatch = true; 
            }
            else if constexpr( std::is_convertible<decltype(arg), includePythonBindings>::value ){
              m_includePythonBindings = true; 
            }
            else if constexpr( std::is_convertible<decltype(arg), includeParameters>::value ) {
              m_includeParameters = true; 
            }
            else ERROR("Unrecognised argument: " << type_string(arg) ); 
          }; 
          for_each( std::tuple<const namedArgs&...>(args...), process_argument);
          if( mps == nullptr ){ DEBUG("No minuit parameterset linked."); }
          resolve(mps);
          if constexpr(std::is_same<ret_type,void>::value ) 
          {
            typedef typename std::remove_pointer<zeroType<arg_types...>>::type zt; 
            m_outputSize = detail::size_of<zt>::value;
            DEBUG( "one element: " << m_outputSize << type_string<zt>() );
          }
          else if constexpr( isVector<ret_type>::value ){
            m_outputSize = detail::size_of<typename ret_type::value_type>::value;    
          }
          else {
            m_outputSize = detail::size_of<ret_type>::value; 
          }
          if( is<TensorExpression>(expression) ){
            m_outputSize *= cast<TensorExpression>(expression).tensor().size();
          }
        }
        
        CompiledExpression( const std::string& name = "" ) : CompiledExpressionBase( name ) { m_outputSize = detail::size_of<ret_type>::value; };
        
        void setDebug( const DebugSymbols& db ){ m_db = db; } 
        
        std::vector<real_t> externBuffer() const override { return m_externals ; } 
        std::string returnTypename() const override { return type_string<ret_type>(); }
        bool use_rto() const override { return std::is_same<ret_type, void>::value; }
        std::vector<std::string> types() const override { return typelist<arg_types...>();}
        void resolve( const MinuitParameterSet* mps=nullptr ){ CompiledExpressionBase::resolve(mps); }
        void setExternals( const std::vector<double>& external ) { m_externals = external; }

        unsigned int getNParams() const { return m_externals.size(); }

        void print() const override 
        {
          INFO( "Name     = " << name() );
          INFO( "Hash     = " << hash() );
          INFO( "IsReady? = " << isReady() << " IsLinked? " << (m_fcn.isLinked() ) );
          INFO( "args     = ["<< vectorToString( m_externals, ", ") <<"]");
          auto func = orderedCacheFunctors();
          for( auto& c : func ) c->print() ;
        }

        void setExternal( const double& value, const unsigned int& address ) override
        {
          if ( m_externals[address] == value ) return;
          DEBUG( "Setting external " << address << " / " << m_externals.size() << " to value = " << value << " ; current = " << m_externals[address] );
          m_externals[address]  = value;
          m_hasExternalsChanged = true;
        }
        void resizeExternalCache(const size_t& N ) override 
        { 
          if( m_externals.size() < N ) m_externals.resize(N);
        }
        bool hasExternalsChanged() { return m_hasExternalsChanged; }
        void resetExternals() { m_hasExternalsChanged = false; }

        const Expression& expression() const { return m_obj; }
        bool isReady()          const override { return m_fcn.isLinked(); }
        bool isLinked()         const { return m_fcn.isLinked() ; } 

        unsigned returnTypeSize() const override { return m_outputSize; }

        template < typename T > ret_type operator()( const T* event ) const
        {
          return m_fcn( m_externals.data(), event );
        }
        ret_type operator()( const arg_types&... args ) const 
        {
          return m_fcn( args... );
        }
        template <typename... batch_arg_types> void batch( batch_arg_types... args ) const { 
          m_batchFcn(args...); 
        }

        template < typename T> void debug( const T* event ) const
        {
          if ( !m_fcn.isLinked() ) {
            FATAL( "Function " << name() << " not linked" );
          }
          if ( !m_fdb.isLinked() ) {
            FATAL( "Function" << name() << " debugging symbols not linked" );
          }
          std::vector<std::pair<std::string, complex_v>> debug_results;
          if constexpr(std::is_same<void, ret_type>::value) debug_results = m_fdb( nullptr, 0, &( m_externals[0] ), event );
          else debug_results = m_fdb( &(m_externals[0]), event);
          for( auto& debug_result : debug_results ){ 
            auto val = debug_result.second;  
            auto label = debug_result.first; 
            if( utils::all_of(val.real(), -999.) )  std::cout << bold_on << std::setw(50) << std::left << label << bold_off << std::endl; 
            else if( utils::all_of(val.imag(), 0.) ) std::cout << "  "    << std::setw(50) << std::left << label << " = " << val.real() << std::endl; 
            else
              std::cout << "  "    << std::setw(50) << std::left << label << " = " << val << std::endl; 
          }
        }

        bool link( void* handle ) override
        {
          const std::string symbol = progName() ;
          bool status = true;
          status &= m_fcn.set(handle, symbol, true);
          status &= m_db.size() == 0 || m_fdb.set(handle, symbol + "_DB");
          if( !m_disableBatch ) status &= m_batchFcn.set(handle, symbol + "_batch");
          return status;
        }
        bool link( const std::string& handle ) override
        {
          return link( dlopen( handle.c_str(), RTLD_NOW ) );
        };
        std::string arg_type(const unsigned& i ) const override
        {
          return typelist<arg_types...>()[i];
        }
    };

  template <typename return_type> 
    CompiledExpression<void(return_type*, const double*, const double*)> 
    make_rto_expression( const Expression& expression, const std::string& name)
    {
      CompiledExpression<void(return_type*, const double*, const double*)> rt(expression,name);
      rt.compile();
      rt.prepare();
      return rt;
    }

  template <typename return_type> CompiledExpression<return_type(const double*, const double*)> 
    make_expression( const Expression& expression, const std::string& name)
    {
      CompiledExpression<return_type(const double*, const double*)> rt(expression,name);
      rt.compile();
      rt.prepare();
      return rt;
    }
  template <typename return_type, typename arg1 = double, typename arg2 =double, typename... arg_types> 
    CompiledExpression<return_type(const arg1*, const arg2*)> 
    make_expression( const Expression& expression, 
        const std::string& name, 
        const arg_types&... args)
    {
      CompiledExpression<return_type(const arg1*, const arg2*)> rt(expression,name,args...);
      rt.compile();
      rt.prepare();
      return rt;
    }

} // namespace AmpGen


#endif
