#include "AmpGen/Integrator2.h"

using namespace AmpGen;

void Integrator2::integrateBlock()
{
  real_t re[N] = {0};
  real_t im[N] = {0};
  size_t addr_i[N] = {0};
  size_t addr_j[N] = {0};
  for( size_t roll = 0 ; roll < N; ++roll )
  {
    addr_i[roll] = m_integrals[roll].i;
    addr_j[roll] = m_integrals[roll].j;
  }
  for ( size_t roll = 0; roll < N; ++roll ) {
    auto& b1 = m_buffer[m_integrals[roll].i];
    auto& b2 = m_buffer[m_integrals[roll].j];
    #pragma omp parallel for reduction(+: re, im)
    for ( size_t i = 0; i < m_events->size(); ++i ) {
      auto c = b1[i] * std::conj(b2[i]);
      re[roll] += m_weight[i] * std::real(c);
      im[roll] += m_weight[i] * std::imag(c);
    }
  }
  for ( size_t j = 0; j < m_counter; ++j ) m_integrals[j].transfer( complex_t( re[j], im[j] ) / m_norm );
  m_counter = 0;
}

Integrator2::Integrator2( const EventList* events ) : m_events( events )
{
  if( m_events == nullptr ) return;
  m_weight.resize( m_events->size() );
  for( size_t i = 0 ; i < m_events->size(); ++i )
  {
    m_weight[i] = m_events->at(i).weight() / m_events->at(i).genPdf();
    m_norm      += m_weight[i]; 
  }
}

bool Integrator2::isReady()            const { return m_events != nullptr; }
const EventList& Integrator2::events() const { return *m_events; } 
void Integrator2::queueIntegral(const size_t& c1, 
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
void Integrator2::addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc )
{
  m_integrals[m_counter++] = Integral<arg>(c1,c2,tFunc);
  if ( m_counter == N ) integrateBlock(); 
}
void Integrator2::queueIntegral(const size_t& i, const size_t& j, complex_t* result)
{
  addIntegralKeyed(i, j, [result](arg& val){ *result = val ; } ); 
}
void Integrator2::flush()
{
  if ( m_counter == 0 ) return;
  integrateBlock();
}
void Integrator2::setBuffer( complex_t* pos, const complex_t& value, const size_t& size )
{
  *pos = value;
}

void Integrator2::setBuffer( complex_t* pos, const std::vector<complex_t>& value, const size_t& size)
{
  memcpy( pos, &(value[0]), size * sizeof(complex_t) );
}
