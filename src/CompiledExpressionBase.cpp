#include "AmpGen/CompiledExpressionBase.h"

#include <algorithm>
#include <memory>
#include <ostream>

#include "AmpGen/CacheTransfer.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ThreadPool.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/ASTResolver.h"

using namespace AmpGen;

std::string AmpGen::programatic_name( std::string s )
{
  std::replace( s.begin(), s.end(), '-', 'm' );
  std::replace( s.begin(), s.end(), '+', 'p' );
  std::replace( s.begin(), s.end(), '-', 'm' );
  std::replace( s.begin(), s.end(), ',', '_' );
  std::replace( s.begin(), s.end(), '(', '_' );
  std::replace( s.begin(), s.end(), ')', '_' );
  std::replace( s.begin(), s.end(), '[', '_' );
  std::replace( s.begin(), s.end(), ']', '_' );
  std::replace( s.begin(), s.end(), '{', '_' );
  std::replace( s.begin(), s.end(), '}', '_' );
  std::replace( s.begin(), s.end(), '*', 's' );
  std::replace( s.begin(), s.end(), ':', '_' );
  std::replace( s.begin(), s.end(), '\'', '_' );
  return s;
}

void CompiledExpressionBase::resolve(const MinuitParameterSet* mps)
{
  ASTResolver resolver( m_evtMap, mps );
  resolver.mps = mps;
  resolver.useCompileTimeConstants = 
    NamedParameter<bool>( "CompiledExpression::UseCompileTimeConstants", false ); 
  m_obj.resolve( resolver );
  std::map< uint64_t , unsigned int> used_functions;
  do { 
    resolver.reduceSubTrees();
    auto subexpressions = resolver.subTrees;
    for( auto& expression : subexpressions ) {     
      expression.second.resolve( resolver );
      auto stack_pos = used_functions.find( expression.first );
      if ( stack_pos == used_functions.end() ) {
        m_dependentSubexpressions.emplace_back( "v"+std::to_string( expression.first ) , expression.second );
        used_functions[expression.first] = m_dependentSubexpressions.size() - 1;
        continue;
      }
      unsigned int oldPos = stack_pos->second;
      auto it             = m_dependentSubexpressions.begin() + oldPos;
      if ( it == m_dependentSubexpressions.end() - 1 ) continue;

      std::rotate( it, it + 1, m_dependentSubexpressions.end() );

      for ( auto uf = used_functions.begin(); uf != used_functions.end(); ++uf ) {
        if ( uf->second >= oldPos ) uf->second = uf->second - 1;
      }
      used_functions[expression.first] = m_dependentSubexpressions.size() - 1;
    }
  } while ( resolver.hasSubExpressions() );
  std::reverse( m_dependentSubexpressions.begin(), m_dependentSubexpressions.end() );

  for ( auto& sym : m_db ) {
    sym.second.resolve( resolver );
    sym.second.resolve( resolver );
  }
  resolveParameters(resolver);
}

CompiledExpressionBase::CompiledExpressionBase( const Expression& expression, 
    const std::string& name, 
    const DebugSymbols& db,
    const std::map<std::string, size_t>& evtMapping )
  : m_obj( expression ), 
  m_name( name ),
  m_db(db),
  m_evtMap(evtMapping),
  m_readyFlag(nullptr) {}

CompiledExpressionBase::CompiledExpressionBase( const std::string& name ) 
  : m_name( name ),
  m_readyFlag(nullptr) {}

  std::string CompiledExpressionBase::name() const { return m_name; }
  unsigned int CompiledExpressionBase::hash() const { return FNV1a_hash(m_name); }

void CompiledExpressionBase::resolveParameters( ASTResolver& resolver )
{
  auto ppdf = this;
  m_cacheTransfers.clear();
  for( auto& expression : resolver.cacheFunctions ) 
    m_cacheTransfers.emplace_back( expression.second ); 
  ppdf->resizeExternalCache(resolver.nParameters); 
  prepare();
}

void CompiledExpressionBase::prepare()
{
  for ( auto& t : m_cacheTransfers ) t->transfer( this );
}

void CompiledExpressionBase::addDependentExpressions( std::ostream& stream, size_t& sizeOfStream ) const
{
  for ( auto& dep : m_dependentSubexpressions ) {
    std::string rt = "auto " + dep.first + " = " + dep.second.to_string() +";"; 
    stream << rt << "\n";
    sizeOfStream += sizeof(char) * rt.size(); /// bytes /// 
  }
}

void CompiledExpressionBase::to_stream( std::ostream& stream  ) const 
{
  if( m_db.size() !=0 ) stream << "#include<iostream>\n"; 
  stream << "extern \"C\" const char* " << m_name << "_name() {\n"
    << "  return \"" << m_name << "\";\n";
  stream << "}\n\n";

  // Avoid a warning about std::complex not being C compatible (it is)
  stream << "#pragma clang diagnostic push\n"
    << "#pragma clang diagnostic ignored \"-Wreturn-type-c-linkage\"\n";

  stream << "extern \"C\" " << returnTypename() << " " << m_name << "(" << fcnSignature() << "){\n";
  size_t sizeOfStream = 0; 
  addDependentExpressions( stream , sizeOfStream );
  std::string objString = m_obj.to_string();
  stream << "return " << objString << ";\n}\n";

  if( NamedParameter<bool>("CompiledExpressionBase::Compat") == true ){
    stream << "#pragma clang diagnostic pop\n\n";
    stream << "extern \"C\" void " <<  m_name << "_c" << "(double *real, double *imag, " << fcnSignature() << "){\n";
    stream << "  auto val = " << m_name << "(" << args() << ") ;\n"; 
    stream << "  *real = val.real();\n";
    stream << "  *imag = val.imag();\n";
    stream << "}\n";
  }
  if ( m_db.size() != 0 ) addDebug( stream );
}

std::ostream& AmpGen::operator<<( std::ostream& os, const CompiledExpressionBase& expression )
{
  expression.to_stream(os);
  return os; 
}

void CompiledExpressionBase::compile(const std::string& fname, const bool& wait  )
{
  if( ! wait ){
  m_readyFlag = new std::shared_future<bool>( ThreadPool::schedule([this,fname](){
        CompilerWrapper().compile(*this,fname);
        return true;} ) );
  }
  else {
    CompilerWrapper(false).compile(*this,fname );
  }
}

void CompiledExpressionBase::addDebug( std::ostream& stream ) const
{
  stream << "extern \"C\" void " << m_name << "_DB(" << fcnSignature() << "){\n";
  size_t sizeOf = 0 ; 
  addDependentExpressions( stream, sizeOf );
  for ( auto& f : m_db ) {
    std::string name      = f.first;
    const auto expression = f.second;
    if ( expression.to_string() == "NULL" )
      stream << "std::cout << \"\\033[1m\" << \"" << name << "\" << \"\\033[0m\" << std::endl ;" << std::endl;
    else {
      name.resize( 20, ' ' );
      stream << "std::cout << \"  \" << \"" << name << "\" << \"=  \"<< " << expression.to_string() << "<< std::endl ;"
        << std::endl;
    }
  }
  stream << "}" << std::endl;
}
