#include "AmpGen/MatrixElement.h"
#if ENABLE_AVX 
  #include "AmpGen/EventListSIMD.h"
#else
  #include "AmpGen/EventList.h"
#endif

using namespace AmpGen; 

MatrixElement::MatrixElement(const Particle& dt, 
    const TotalCoupling& coupling, 
    const amp_type& amp) : 
  amp_type(amp),
  decayTree(dt), 
  coupling(coupling) {
        
    if( is<TensorExpression>(amp.expression() ) )
    {
      size = cast<TensorExpression>( amp.expression() ) . size() ; 
    }
  }

MatrixElement::MatrixElement(const Particle& dt, 
    const TotalCoupling& coupling, 
    const MinuitParameterSet& mps,
    const std::map<std::string, unsigned>& evtFormat, 
    const bool& debugThis) :
  amp_type(dt.decayDescriptor()),
  decayTree(dt),
  coupling(coupling)
{
  this->set(Particle(dt).getExpression(debugThis ? &db : nullptr ), dt.decayDescriptor(), evtFormat, db, &mps );
}

std::vector<size_t> AmpGen::processIndex(const std::vector<MatrixElement>& tm, const std::string& label)
{
  std::vector<size_t> indices;
  for ( size_t i = 0; i < tm.size(); ++i ) {
    if ( tm[i].coupling.contains(label) ) indices.push_back(i);
  }
  return indices;
}

const std::vector<complex_v> MatrixElement::operator()(const Event& event) const 
{ 
  std::vector<complex_v> rt(size); 
  std::vector<size_t> offsets(size);
  for( int i = 0 ; i != size; ++i ) offsets[i] = i;
  amp_type::operator()(rt.data(), offsets.data(), externBuffer().data(), 
      #if ENABLE_AVX 
      EventListSIMD::makeEvent(event).data()
      #else
      event.address() 
      #endif 
      ); 
  return rt;
}

size_t AmpGen::findIndex(const std::vector<MatrixElement>& tm, const std::string& decayDescriptor)
{
  for ( size_t i = 0; i < tm.size(); ++i ) {
    if ( tm[i].decayDescriptor() == decayDescriptor ) return i;
  }
  ERROR( "Component " << decayDescriptor << " not found" );
  return 999;
}

void MatrixElement::debug( const Event& event ) const 
{
#if ENABLE_AVX
  amp_type::debug(EventListSIMD::makeEvent(event).data() );
#else
  amp_type::debug(event.address()); 
#endif 
}

std::vector<size_t> AmpGen::findIndices(const std::vector<MatrixElement>& tm, const std::string& decayDescriptor)
{
  std::vector<size_t> rt; 
  for ( size_t i = 0; i < tm.size(); ++i ) {
    if ( tm[i].decayDescriptor().find(decayDescriptor) != std::string::npos ) rt.push_back(i);
  }
  return rt;
}
