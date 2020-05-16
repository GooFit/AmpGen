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
#include <cxxabi.h>
#include <dlfcn.h>
#include <vector>

namespace AmpGen
{
  /* @class CompiledExpression
     @tparam RETURN_TYPE The type that is returned this compiled expression,
     usually this is a std::complex<double>,
     but in principal support also exists for computing coupled channel propagators
     (i.e. returning array types) */
  namespace detail {
    template <typename T> struct size_of {       static constexpr unsigned value = sizeof(T); };
    template <>           struct size_of<void> { static constexpr unsigned value = 0; } ;
  }
  template <class RETURN_TYPE, class... ARGS> class CompiledExpression; 

  template <class RETURN_TYPE, class... ARGS>
    class CompiledExpression<RETURN_TYPE(ARGS...)> : public CompiledExpressionBase
    {

      private:
        DynamicFCN<RETURN_TYPE( ARGS... )>                                       m_fcn;
        DynamicFCN<void( const size_t&, const size_t&, const size_t&, RETURN_TYPE*, ARGS... )> m_batchFcn;
        DynamicFCN<std::vector<std::pair<std::string, complex_v>>(ARGS...)>      m_fdb;
        std::vector<real_t>  m_externals             = {};
        bool                 m_hasExternalsChanged   = {false};

      public:
        typedef RETURN_TYPE return_type;

        CompiledExpression( const Expression& expression, 
            const std::string& name,
            const std::map<std::string, unsigned>& evtMapping = 
            std::map<std::string, unsigned>(),
            const DebugSymbols& db = {},
            const MinuitParameterSet* mps = nullptr )
          : CompiledExpressionBase( expression, name, db, evtMapping ) 
        {
          resolve(mps);
        }

        CompiledExpression( const std::string& name = "" ) : CompiledExpressionBase( name ) {};
        std::vector<real_t> externBuffer() const { return m_externals ; } 
        std::string returnTypename() const override { return type_string<RETURN_TYPE>(); }
        std::string fcnSignature() const override
        {
          return CompiledExpressionBase::fcnSignature(typelist<ARGS...>(), use_rto());
        }
        bool use_rto() const override {
          return std::is_same<RETURN_TYPE, void>::value;   
        }
        std::string args() const override 
        {
          std::string signature; 
          auto argTypes = typelist<ARGS...>();
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
          std::vector<const CacheTransfer*> ordered_cache_functors; 
          for( const auto& c : m_cacheTransfers ) ordered_cache_functors.push_back( c.get() );
          std::sort( ordered_cache_functors.begin(),
                     ordered_cache_functors.end(),
                     [](auto& c1, auto& c2 ) { return c1->address() < c2->address() ; } );
          for( auto& c : ordered_cache_functors ) c->print() ;
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
          stream << "#include <omp.h>\n";
          stream << "extern \"C\" void " << progName() 
                 << "_batch(";
          stream << " const size_t& N, " 
                 << " const size_t& eventSize, " 
                 << " const size_t& cacheSize, ";
          stream <<  type_string<return_type>() << " * rt, ";
          stream << CompiledExpressionBase::fcnSignature(typelist<ARGS...>(), use_rto(), false) << ") {\n";
          stream << "#pragma omp parallel for\n";
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
          stream << "  double externalParameters [] = {" << (m_externals.size() == 0 ? "0" : vectorToString(m_externals,", ") ) <<"};\n" ;
          stream << "  return " << progName() << "( externalParameters, E ); // E is P \n}\n";
        }

        bool isReady()          const override { return m_fcn.isLinked(); }
        bool isLinked()         const { return m_fcn.isLinked() ; } 

        unsigned returnTypeSize() const override { return detail::size_of<RETURN_TYPE>::value; }

        template < class T > 
          RETURN_TYPE operator()( const T* event ) const
          {
            return m_fcn( m_externals.data(), event );
          }
        RETURN_TYPE operator()( const ARGS&... args ) const 
        {
          return m_fcn( args... );
        }
        template <class... arg_types> void batch( arg_types... args ) const { 
          m_batchFcn(args...); 
        }

        template < class T> void debug( const T* event ) const
        {
          if ( !m_fcn.isLinked() ) {
            FATAL( "Function " << name() << " not linked" );
          }
          if ( !m_fdb.isLinked() ) {
            FATAL( "Function" << name() << " debugging symbols not linked" );
          }
          std::vector<std::pair<std::string, complex_v>> debug_results;
          if constexpr(std::is_same<void, RETURN_TYPE>::value) debug_results = m_fdb( nullptr, 0, &( m_externals[0] ), event );
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
          status &= !m_enableBatch   || m_batchFcn.set(handle, symbol + "_batch");
          return status;
        }
        bool link( const std::string& handle ) override
        {
          return link( dlopen( handle.c_str(), RTLD_NOW ) );
        };
    };

  template <class RT> 
    CompiledExpression<void(RT*, const double*, const double*)> 
    make_rto_expression( const Expression& expression, const std::string& name , const bool& verbose=false)
    {
      CompiledExpression<void(RT*, const double*, const double*)> rt(expression,name);
      rt.compile();
      rt.prepare();
      return rt;
    }

  template <class RT> 
    CompiledExpression<RT(const double*, const double*)> 
    make_expression( const Expression& expression, const std::string& name , const bool& verbose=false)
    {
      CompiledExpression<RT(const double*, const double*)> rt(expression,name);
      rt.compile();
      rt.prepare();
      return rt;
    }
  template <class RT, class arg1 = double, class arg2 =double> 
    CompiledExpression<RT(const arg1*, const arg2*)> 
    make_expression( const Expression& expression, 
        const std::string& name, 
        const MinuitParameterSet& mps )
    {
      CompiledExpression<RT(const arg1*, const arg2*)> rt(expression,name,{},{},&mps);
      rt.compile();
      rt.prepare();
      return rt;
    }
  template <class RT> 
    CompiledExpression<RT(const double*, const double*)> 
    make_expression( const Expression& expression, 
        const std::string& name,
        const std::map<std::string, unsigned> & evtMap,
        const MinuitParameterSet& mps )
    {
      CompiledExpression<RT(const double*, const double*)> rt(expression,name,evtMap,{},&mps);
      rt.compile();
      rt.prepare();
      return rt;
    }
} // namespace AmpGen


#endif
