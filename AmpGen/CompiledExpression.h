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

        template <typename... namedArgs> 
        CompiledExpression( const Expression& expression, const std::string& name, const namedArgs&... args ) : CompiledExpressionBase(expression, name) 
        {
          const MinuitParameterSet* mps = nullptr; 
          auto process_argument = [this, &mps]( const auto& arg ) mutable
          { 
            if constexpr( std::is_convertible<decltype(arg), DebugSymbols>::value  ) this->m_db = arg;
            else if constexpr( std::is_convertible<decltype(arg), std::map<std::string, unsigned>>::value ) this->m_evtMap = arg;
            else if constexpr( std::is_convertible<decltype(arg), const MinuitParameterSet*>::value or 
                               std::is_convertible<decltype(arg), const AmpGen::MinuitParameterSet*>::value or
                               std::is_convertible<decltype(arg), MinuitParameterSet*>::value ){
              mps = arg;
            }
            else if constexpr( std::is_convertible<decltype(arg), MinuitParameterSet>::value ) mps = &arg;
            
            else if constexpr( std::is_convertible<decltype(arg), disableBatch>::value ) {
              WARNING("Disabling bulk evaluation: did you do this on purpose?");
              m_disableBatch = true; 
            }
            else ERROR("Unrecognised argument: " << type_string(arg) ); 
          }; 
          for_each( std::tuple<const namedArgs&...>(args...), process_argument);
          if( mps == nullptr ) WARNING("No minuit parameterset linked.");
          resolve(mps);
          if constexpr(std::is_same<ret_type,void>::value ) 
          {
            typedef typename std::remove_pointer<zeroType<arg_types...>>::type zt; 
            m_outputSize = detail::size_of<zt>::value;
            DEBUG( "one element: " << m_outputSize << type_string<zt>() );
            if( is<TensorExpression>(expression) ) m_outputSize *= cast<TensorExpression>(expression).tensor().size();
          }
          else m_outputSize = detail::size_of<ret_type>::value; 
        }

        CompiledExpression( const std::string& name = "" ) : CompiledExpressionBase( name ) { m_outputSize = detail::size_of<ret_type>::value; };
        std::vector<real_t> externBuffer() const { return m_externals ; } 
        std::string returnTypename() const override { return type_string<ret_type>(); }
        std::string fcnSignature() const override
        {
          return CompiledExpressionBase::fcnSignature(typelist<arg_types...>(), use_rto());
        }
        bool use_rto() const override {
          return std::is_same<ret_type, void>::value;   
        }
        std::string args() const override 
        {
          std::string signature; 
          auto argTypes = typelist<arg_types...>();
          for( unsigned int i = 0 ; i < argTypes.size(); ++i )
          {
            signature += " x"+std::to_string(i) ;
            if( i != argTypes.size() - 1 ) signature += ", ";
          }
          return signature;
        }

        void resolve( const MinuitParameterSet* mps=nullptr )
        {
          CompiledExpressionBase::resolve(mps);
        }
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
          DEBUG( "Setting external " << address << " / " << m_externals.size() << " to value = " << value << " ; current = " << m_externals[address] );
          if ( m_externals[address] == value ) return;
          m_externals[address]  = value;
          m_hasExternalsChanged = true;
        }
        void resizeExternalCache(const size_t& N ) override 
        { 
          if( m_externals.size() < N ) m_externals.resize(N);
        }
        bool hasExternalsChanged() { return m_hasExternalsChanged; }
        void resetExternals() { m_hasExternalsChanged = false; }

        Expression& expression() { return m_obj; }

        void compileDetails( std::ostream& stream ) const
        {
          stream << "extern \"C\" int " << progName() << "_pSize () {\n"
            << "  return " << m_externals.size() << ";\n";
          stream << "}\n";

          stream << "extern \"C\" double " << progName() << "_pVal (int n) {\n";
          for ( size_t i = 0; i < m_externals.size(); i++ )
            stream << "  if(n == " << i << ") return  " << m_externals.at( i ) << ";\n";
          stream << "  return 0;\n}\n";
        }
        void compileBatch( std::ostream& stream ) const override  
        {
          #if USE_OPENMP
          stream << "#include <omp.h>\n";
          #endif
          stream << "extern \"C\" void " << progName() 
                 << "_batch(";
          stream << " const size_t& N, " 
                 << " const size_t& eventSize, " 
                 << " const size_t& cacheSize, ";
          stream <<  type_string<ret_type>() << " * rt, ";
          stream << CompiledExpressionBase::fcnSignature(typelist<arg_types...>(), use_rto(), false) << ") {\n";
          #if USE_OPENMP
          stream << "#pragma omp parallel for\n";
          #endif
          stream << "for( size_t i = 0; i < N/" << utils::size<float_v>::value << "; ++i ){\n";
          if( use_rto() ) stream << progName() + "( r + cacheSize * i, s, x0, x1 +  i * eventSize);";
          else            stream << " rt[cacheSize*i] = " << progName() + "( x0, x1 +  i * eventSize);";
          stream << "}\n}";
        }

        void compileWithParameters( std::ostream& stream ) const
        {
          DEBUG( "Compiling " << name() << " = " << hash() );
          stream << "extern \"C\" " << returnTypename() << " " << progName() << "_wParams"
            << "( const double*__restrict__ E ){" << std::endl;
          stream << "  double externalParameters [] = {" << 
            (m_externals.size() == 0 ? "0" : vectorToString(m_externals,", ", [](auto& line){ return std::to_string(line); }) ) <<"};\n" ;
          stream << "  return " << progName() << "( externalParameters, E ); // E is P \n}\n";
        }

        bool isReady()          const override { return m_fcn.isLinked(); }
        bool isLinked()         const { return m_fcn.isLinked() ; } 

        unsigned returnTypeSize() const override { return m_outputSize; }

        template < typename T > 
          ret_type operator()( const T* event ) const
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
    make_rto_expression( const Expression& expression, const std::string& name , const bool& verbose=false)
    {
      CompiledExpression<void(return_type*, const double*, const double*)> rt(expression,name);
      rt.compile();
      rt.prepare();
      return rt;
    }

  template <typename return_type> CompiledExpression<return_type(const double*, const double*)> 
    make_expression( const Expression& expression, const std::string& name , const bool& verbose=false)
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
