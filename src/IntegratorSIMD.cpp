#include "AmpGen/IntegratorSIMD.h"
#include "AmpGen/simd/utils.h"

using namespace AmpGen;
using namespace AmpGen::AVX2;

void IntegratorSIMD::integrateBlock()
{
  #pragma omp parallel for
  for ( size_t roll = 0; roll < N; ++roll ) {
    float_v re( _mm256_set1_ps(0.f) );
    float_v im( _mm256_set1_ps(0.f) );
    auto b1 = m_cache.data() + m_integrals[roll].i * m_events->size();
    auto b2 = m_cache.data() + m_integrals[roll].j * m_events->size();
    for ( size_t i = 0; i < m_events->nBlocks(); ++i ) {
      auto c = b1[i] * conj(b2[i]);
      re = _mm256_fmadd_ps(re, m_weight[i], real(c) );
      im = _mm256_fmadd_ps(im, m_weight[i], imag(c) );
    }
    m_integrals[roll].transfer( complex_t( utils::sum_elements(float_v(re)), 
                                           utils::sum_elements(float_v(im)) ) / m_norm );
  }
  m_counter = 0;
}

IntegratorSIMD::IntegratorSIMD( const EventListSIMD* events ) : m_events( events )
{
  if( m_events == nullptr ) return;
  m_weight.resize( m_events->nBlocks() );
  float_v norm_acc = 0.;
  for( size_t i = 0 ; i < m_events->nBlocks(); ++i )
  {
    m_weight[i]  = m_events->weight(i) / m_events->genPDF(i);
    norm_acc = norm_acc + m_weight[i]; 
  }
  m_norm = utils::sum_elements(norm_acc);
}

bool IntegratorSIMD::isReady()            const { return m_events != nullptr; }

void IntegratorSIMD::queueIntegral(const size_t& c1, 
    const size_t& c2, 
    const size_t& i, 
    const size_t& j, 
    Bilinears* out, 
    const bool& sim)
{
  if( !out->workToDo(i,j) ) return;
  if( sim ) 
    addIntegralKeyed( c1, c2, [out,i,j]( arg& val ){ 
        out->set(i,j,val);
        if( i != j ) out->set(j,i, std::conj(val) ); } );
  else 
    addIntegralKeyed( c1, c2, [out,i,j]( arg& val ){ out->set(i,j,val); } );
}
void IntegratorSIMD::addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc )
{
  m_integrals[m_counter++] = Integral<arg>(c1,c2,tFunc);
  if ( m_counter == N ) integrateBlock(); 
}
void IntegratorSIMD::queueIntegral(const size_t& i, const size_t& j, complex_t* result)
{
  addIntegralKeyed(i, j, [result](arg& val){ *result = val ; } ); 
}
void IntegratorSIMD::flush()
{
  if ( m_counter == 0 ) return;
  integrateBlock();
}
