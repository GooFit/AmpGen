#include "AmpGen/CompiledExpressionBase.h"

#include <cxxabi.h>
#include <algorithm>
#include <memory>
#include <ostream>

#include "AmpGen/CacheTransfer.h"
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
  std::replace( s.begin(), s.end(), '*', 's' );
  std::replace_if( s.begin(), s.end(), [](auto& c){ return ! std::isalnum(c) ; }, '_' );
  if( isdigit( s[0] ) ) s = "f" + s;
  std::replace( s.begin(), s.end(), '\'', '_' );
  return s;
}

void CompiledExpressionBase::resolve(const MinuitParameterSet* mps)
{
  if( m_resolver != nullptr ) delete m_resolver ; 
  m_resolver = new ASTResolver( m_evtMap, mps );
  m_resolver->getOrderedSubExpressions( m_obj,  m_dependentSubexpressions ); 
  for ( auto& sym : m_db ){ 
    sym.second.resolve( *m_resolver ); 
    m_resolver->getOrderedSubExpressions( sym.second, m_debugSubexpressions );
  }
  m_cacheTransfers.clear();
  for( auto& expression : m_resolver->cacheFunctions() ) 
    m_cacheTransfers.emplace_back( expression.second ); 
  resizeExternalCache( m_resolver->nParams() ); 
  prepare(); 
}

CompiledExpressionBase::CompiledExpressionBase( const Expression& expression, 
    const std::string& name, 
    const DebugSymbols& db,
    const std::map<std::string, size_t>& evtMapping )
  : m_obj( expression ), 
  m_name( name ),
  m_progName( programatic_name(name) ),
  m_db(db),
  m_evtMap(evtMapping),
  m_readyFlag(nullptr) {}

CompiledExpressionBase::CompiledExpressionBase( const std::string& name ) 
  : m_name( name ),
    m_progName( programatic_name(name) ),
    m_readyFlag(nullptr) {}

std::string CompiledExpressionBase::name() const { return m_name; }
std::string CompiledExpressionBase::progName() const { return m_progName; }
unsigned int CompiledExpressionBase::hash() const { return FNV1a_hash(m_name); }

void CompiledExpressionBase::prepare()
{
  for ( auto& t : m_cacheTransfers ) t->transfer( this );
}

void CompiledExpressionBase::addDependentExpressions( std::ostream& stream, size_t& sizeOfStream ) const
{
  for ( auto& dep : m_dependentSubexpressions ) {
    std::string rt = "auto v" + std::to_string(dep.first) + " = " + dep.second.to_string(m_resolver) +";"; 
    stream << rt << "\n";
    sizeOfStream += sizeof(char) * rt.size(); /// bytes /// 
  }
}

void CompiledExpressionBase::to_stream( std::ostream& stream  ) const 
{
  if( m_db.size() !=0 ) stream << "#include<iostream>\n"; 
  stream << "extern \"C\" const char* " << progName() << "_name() {  return \"" << m_name << "\"; } \n";
  bool enable_cuda = NamedParameter<bool>("EnableCUDA",false);

    size_t sizeOfStream = 0; 
  if( !enable_cuda ){
    // Avoid a warning about std::complex not being C compatible (it is)
    stream << "#pragma clang diagnostic push\n"
      << "#pragma clang diagnostic ignored \"-Wreturn-type-c-linkage\"\n";

    stream << "extern \"C\" " << returnTypename() << " " << progName() << "(" << fcnSignature() << "){\n";
    addDependentExpressions( stream , sizeOfStream );
    std::string objString = m_obj.to_string(m_resolver);
    stream << "return " << objString << ";\n}\n";
  }
  else {
    std::string rt_cpp = returnTypename();
    std::string rt_cuda = "";
    if( rt_cpp == "double"               || rt_cpp == "float"               || rt_cpp == "real_t" ) 
      rt_cuda = "float* r, const int N";
    if( rt_cpp == "std::complex<double>" || rt_cpp == "std::complex<float>" || rt_cpp == "complex_t" ) 
      rt_cuda = "ampgen_cuda::complex_t* r, const int N";
    stream << "__global__ void " << progName() << "( " << rt_cuda << ", const float_t* x0, const float3* x1){\n"; 
    stream <<  "  int i     = blockIdx.x * blockDim.x + threadIdx.x;\n";
    addDependentExpressions( stream, sizeOfStream);
    std::string objString = m_obj.to_string(m_resolver);
    stream << "  r[i] = " << objString << ";\n}\n";
  }

  if( NamedParameter<bool>("IncludePythonBindings", false) == true ){
    stream << "#pragma clang diagnostic pop\n\n";
    stream << "extern \"C\" void " <<  progName() << "_c" << "(double *real, double *imag, " << fcnSignature() << "){\n";
    stream << "  auto val = " << progName() << "(" << args() << ") ;\n"; 
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
  if(!wait){
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
  stream << "#include<string>\n";
  stream << "extern \"C\" std::vector<std::pair< std::string, std::complex<double>>> " 
         << m_progName << "_DB(" << fcnSignature() << "){\n";
  for ( auto& dep : m_debugSubexpressions ) {
    std::string rt = "auto v" + std::to_string(dep.first) + " = " + dep.second.to_string(m_resolver) +";"; 
    stream << rt << "\n";
  }
  stream << "return {";
  for ( unsigned int i = 0 ; i < m_db.size(); ++i ) {
    std::string comma = (i!=m_db.size()-1)?", " :"};\n}\n";
    const auto expression = m_db[i].second; 
    stream << std::endl << "{\"" << m_db[i].first << "\",";
    if ( expression.to_string(m_resolver) != "NULL" )
      stream << expression.to_string(m_resolver) << "}" << comma;
    else stream << "-999}" << comma ;
  }
}
