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

  template <class RETURN_TYPE, class... ARGS>
  class CompiledExpression : public CompiledExpressionBase
  {

  private:
    DynamicFCN<RETURN_TYPE( ARGS... )>                                  m_fcn;
    DynamicFCN<std::vector<std::pair<std::string, complex_t>>(ARGS...)> m_fdb;
    std::vector<real_t>  m_externals             = {};
    bool                 m_hasExternalsChanged   = {false};
   
  public:

    CompiledExpression( const Expression& expression, 
                        const std::string& name,
                        const std::map<std::string, size_t>& evtMapping = 
                              std::map<std::string, size_t>(),
                        const DebugSymbols& db = {},
                        const MinuitParameterSet* mps = nullptr )
      : CompiledExpressionBase( expression, name, db, evtMapping ) 
    {
      resolve(mps);
    }

    CompiledExpression( const std::string& name = "" ) : CompiledExpressionBase( name ) {};
    std::vector<real_t> externBuffer() const { return m_externals ; } 
    std::string returnTypename() const override { return typeof<RETURN_TYPE>(); }
    std::string fcnSignature() const override
    {
      std::string signature;
      auto argTypes = typelist<ARGS...>();
      for( unsigned int i = 0 ; i < argTypes.size(); ++i )
      {
        signature += argTypes[i] + " x"+std::to_string(i) ;
        if( i != argTypes.size() - 1 ) signature += ", ";
      }
      return signature;
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
      for( auto& c : m_cacheTransfers ){ c->print() ; } 
    }

    void setExternal( const double& value, const unsigned int& address ) override
    {
      DEBUG( "Setting external " << address << " / " << m_externals.size() << " to value = " << value << " ; current = " << m_externals[address] );
      if ( m_externals[address] == value ) return;
      m_externals[address]  = value;
      m_hasExternalsChanged = true;
    }
    void resizeExternalCache(const size_t& N ) override { 
      if( m_externals.size() < N ){
        m_externals.resize(N);
      }
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

    void compileWithParameters( std::ostream& stream ) const
    {
      DEBUG( "Compiling " << name() << " = " << hash() );
      // Avoid a warning about std::complex not being C compatible (it is)
      // stream << "#pragma clang diagnostic push\n"
      //        << "#pragma clang diagnostic ignored \"-Wreturn-type-c-linkage\"\n";
      stream << "extern \"C\" " << returnTypename() << " " << progName() << "_wParams"
             << "( const double*__restrict__ E ){" << std::endl;
      stream << "  double externalParameters [] = {";
      if ( m_externals.size() != 0 ) {
        for ( unsigned int i = 0; i < m_externals.size() - 1; ++i ) {
          stream << m_externals[i] << ", " ; 
        }
        stream << m_externals[m_externals.size() - 1] << " }; " << std::endl;
      } else stream << "0};" << std::endl;
      stream << "  return " << progName() << "( externalParameters, E ); // E is P \n}\n";
      // stream << "#pragma clang diagnostic pop\n\n" << std::endl;
    }

    bool isReady()          const override { return m_fcn.isLinked(); }
    bool isLinked()         const { return m_fcn.isLinked() ; } 
    size_t returnTypeSize() const override { return sizeof( RETURN_TYPE ); }
   
    template < class T > 
    RETURN_TYPE operator()( const T* event ) const
    {
      return m_fcn( (const double*)( &( m_externals[0] ) ), event );
    }
    RETURN_TYPE operator()( const ARGS&... args ) const 
    {
      return m_fcn( args... );
    } 
    template < class T> 
    void debug( const T* event ) const
    {
      if ( !m_fcn.isLinked() ) {
        FATAL( "Function " << name() << " not linked" );
      }
      if ( !m_fdb.isLinked() ) {
        FATAL( "Function" << name() << " debugging symbols not linked" );
      }
      auto debug_results = m_fdb( &( m_externals[0] ), event );
      for( auto& debug_result : debug_results ){ 
        auto val = debug_result.second;  
        auto label = debug_result.first; 
        if( std::real(val) == -999. )  std::cout << bold_on << std::setw(50) << std::left << label << bold_off << std::endl; 
        else if( std::imag(val) == 0 ) std::cout << "  " << std::setw(50) << std::left << label << " = " << std::real(val) << std::endl; 
        else                           std::cout << "  " << std::setw(50) << std::left << label << " = " << val << std::endl; 
      }
    }

    bool link( void* handle ) override
    {
      const std::string symbol = progName() ;
      if ( m_fcn.set( handle, symbol ) == 0 ) {
        ERROR( dlerror() );
        FATAL( name() << " (symbol = " << symbol << ") linking fails" );
        return false;
      }
      if ( m_db.size() ==0 ) return true;
      if ( m_fdb.set( handle, progName() + "_DB" ) == 0 ) {
        ERROR( "Linking of " << name() << " symbol = " << symbol << ") for debugging fails" );
        return false;
      }
      return true;
    }
    bool link( const std::string& handle ) override
    {
      DEBUG( "Linking " << name() << ( m_db.size() !=0  ? " (debugging)" : "" ) << "    hash = " << hash()  );
      const std::string symbol = progName();
      if ( m_fcn.set( handle, symbol ) == 0 ) {
        ERROR( "Function not linked: " << name() << " (sym="<< symbol << ")" );
        return false;
      }
      if ( m_db.size() ==0 ) return true;
      const std::string dbsymbol = symbol + "_DB";
      if ( m_fdb.set( handle, dbsymbol ) == 0 ) {
        ERROR( "Linking of " << name() << " symbol = " << dbsymbol << ")" );
        return false;
      }
      return true;
    }
  };

  template <class RT> 
    CompiledExpression<RT, const double*, const double*> 
      make_expression( const Expression& expression, const std::string& name , const bool& verbose=false)
      {
        CompiledExpression<RT,const double*, const double*> rt(expression,name);
        rt.compile();
        rt.prepare();
        return rt;
      }
  template <class RT> 
    CompiledExpression<RT, const double*, const double*> 
      make_expression( const Expression& expression, 
                       const std::string& name, 
                       const MinuitParameterSet& mps )
      {
        CompiledExpression<RT,const double*, const double*> rt(expression,name,{},{},&mps);
        rt.compile();
        rt.prepare();
        return rt;
      }
  template <class RT> 
    CompiledExpression<RT, const double*, const double*> 
      make_expression( const Expression& expression, 
                       const std::string& name,
                       const std::map<std::string, size_t> & evtMap,
                       const MinuitParameterSet& mps )
      {
        CompiledExpression<RT,const double*, const double*> rt(expression,name,evtMap,{},&mps);
        rt.compile();
        rt.prepare();
        return rt;
      }
} // namespace AmpGen


#endif
