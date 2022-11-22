#include "AmpGen/IncoherentSum.h"

#include <memory.h>
#include <iomanip>
#include <memory>
#include <ostream>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/EventList.h"

using namespace AmpGen;

IncoherentSum::IncoherentSum( const EventType& finalStates, 
                              const MinuitParameterSet& mps,
                              const std::string& prefix ) : 
  CoherentSum( finalStates, mps, prefix )
{
  //m_normalisations.resize( size(), 1 );
}

double IncoherentSum::norm() const
{ 
  return norm( m_normalisations ); 
}

double IncoherentSum::norm( const Bilinears& norms ) const
{
  double norm( 0 ); // (0,0);
  for ( unsigned int i = 0; i < size(); ++i ) {
    double val = norms.get( i, i ).real() * std::norm( m_matrixElements[i].coefficient );
    norm += val;
  }
  return norm; //.real();
}

//void IncoherentSum::prepare()
//{
  /*
  if ( m_isConstant && m_prepareCalls != 0 ) return;
  transferParameters();
  for ( auto& mE : m_matrixElements ) {
    auto& amp = mE.amp;
    amp.prepare();
    if ( m_prepareCalls != 0 && !amp.hasExternalsChanged() ) continue;
    if ( m_prepareCalls == 0 && m_events != nullptr )
      mE.addressData = m_events->registerExpression( amp );
//    if ( m_events != nullptr ) m_events->updateCache( amp, mE.addressData ); 
   // if ( m_prepareCalls == 0 && m_integrator.isReady() ){
   //   m_integrator.prepareExpression( amp );
   // }
    amp.resetExternals();
  }
  if( m_prepareCalls == 0 ){
    for( size_t i = 0 ; i < m_matrixElements.size(); ++i ){
      auto index = m_integrator.events().getCacheIndex( m_matrixElements[i].amp );
      m_integrator.queueIntegral( index, index, i, 0, &m_normalisations, false);
    }
    m_integrator.flush();
  }
  m_prepareCalls++;
  m_norm = norm();
  INFO( "norm = " << m_norm << " weight = " << m_weight );
  */
//}


std::vector<FitFraction> IncoherentSum::fitFractions( const LinearErrorPropagator& linProp )
{
  std::vector<FitFraction> outputFractions;
  std::vector<size_t> normSet(m_matrixElements.size());
  std::iota( std::begin(normSet), std::end(normSet), 0 );
  FitFractionCalculator<IncoherentSum> calc(this, normSet ); 
  for ( unsigned int i = 0; i < m_matrixElements.size(); ++i ) 
    calc.emplace_back(m_matrixElements[i].decayDescriptor(), std::vector<size_t>({i}));

  auto fractions = calc(m_evtType.mother(), linProp);
  for(auto& f : fractions) outputFractions.push_back(f); 

  for ( auto& p : outputFractions ) INFO(p);
  return outputFractions;
}


double IncoherentSum::prob( const Event& evt ) const 
{ 
  return m_weight * prob_unnormalised(evt) / m_norm; 
}

double IncoherentSum::prob_unnormalised( const Event& evt ) const 
{ 
  real_v value( 0. );
  for (unsigned int i = 0 ; i != m_matrixElements.size(); ++i ) {
    value += std::norm( complex_v(m_matrixElements[i].coefficient) * m_cache(evt.index() / utils::size<real_v>::value, i ) );
  }
  #if ENABLE_AVX
    return value.at(evt.index() % utils::size<real_v>::value );
  #else 
    return value;
  #endif      
}


real_v IncoherentSum::operator()( const real_v* /*evt*/, const unsigned block ) const 
{
  real_v value( 0. );
  for ( const auto& mE : m_matrixElements ) 
  {
    unsigned address = &mE - &m_matrixElements[0];
    value += std::norm( complex_v( mE.coefficient) * m_cache(block, address) ); 
  }
  return (m_weight/m_norm ) * value; 
}

#if ENABLE_AVX
double IncoherentSum::operator()( const double* /*evt*/, const unsigned block ) const 
{
  return operator()((const real_v*)nullptr, block / utils::size<real_v>::value ).at( block % utils::size<real_v>::value );
}
#endif


real_t IncoherentSum::prob_unnormalisedNoCache( const Event& evt ) const
{
  return utils::get<0>( real_t(std::accumulate( m_matrixElements.begin(), 
          m_matrixElements.end(), 
          real_t(0), 
          [&evt]( const auto& a, const auto& b ){ auto amp = utils::at(b(evt)[0], 0);
          return a + std::norm(b.coefficient * complex_t(amp));} )) );
}

void IncoherentSum::debug( const Event& evt, const std::string& nameMustContain )
{
  prepare();
  INFO("Weight = " << evt.weight() << " genPDF = " << evt.genPdf() );

  for ( auto& me : m_matrixElements ) {
    auto A = me(evt);
    INFO( std::setw(70) << me.decayTree.uniqueString() 
        << " A = [ "  << A[0].real()             << " " << A[0].imag()
        << " ] g = [ "<< me.coupling().real() << " " << me.coupling().imag() << " ] "
        << m_cache( evt.index(), std::distance(&m_matrixElements[0], &me ) )
        << me.decayTree.CP() );
  }
  if( m_dbThis ) for ( auto& me : m_matrixElements ) me.debug( evt) ;
  INFO( "A(x) = " << prob_unnormalised(evt) << " without cache: " << prob_unnormalisedNoCache(evt) );
}

std::function<real_t(const Event&)> IncoherentSum::evaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  Store<complex_v, Alignment::AoS> store( events->size(), m_matrixElements);
  for( auto& me : m_matrixElements ) store.update(events->store(), me );
  
  std::vector<double> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    real_v amp(0.);
    for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ) 
      amp = amp + std::norm( complex_v(m_matrixElements[j].coefficient) * store(block, j) );
    utils::store( values.data() + block * utils::size<real_v>::value,  (m_weight/m_norm) * amp  );
  }
  return arrayToFunctor<double, typename EventList_type::value_type>(values);
}

KeyedFunctors<double(Event)> IncoherentSum::componentEvaluator(const EventList_type* ievents) const 
{
  using store_t = Store<complex_v, Alignment::SoA>; 
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  KeyedFunctors<double(Event)> rt; 
  std::shared_ptr<const store_t> cache;
  if( events != m_integrator.events<EventList_type>() )
  {
    cache = std::make_shared<const store_t>(events->size(), m_matrixElements);
    for( auto& me : m_matrixElements ) const_cast<store_t*>(cache.get())->update(events->store(), me);
  }
  else cache = std::shared_ptr<const store_t>( & m_integrator.cache(), [](const store_t* t){} ); 
  /// this little slice of weirdness allows either a new cache to be instantiated, or one to just get a pointer to the one used for the integration. 

  for( unsigned i = 0 ; i != m_matrixElements.size(); ++i )
  {
      auto mi = m_matrixElements[i]; 
      auto ci = this->m_matrixElements[i].coefficient;
      auto name = programatic_name(mi.decayTree.decayDescriptor());
      
      auto functor = [ci,i,cache](const Event& event){ return std::norm( ci * cache->get<complex_t>( event.index(), i )  ) ;};
      rt.add(functor, name, "");

  }
  return rt; 
}

