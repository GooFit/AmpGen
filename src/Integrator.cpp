#include "AmpGen/Integrator.h"

using namespace AmpGen;

Bilinears::Bilinears(const size_t& r, const size_t& c) : 
  rows(r), 
  cols(c), 
  norms(c*r), 
  markAsZero(c*r,0), 
  calculate(c*r,1) 
{}

complex_t Bilinears::get( const size_t& x, const size_t& y ) const { return norms[x * cols + y]; }
void      Bilinears::set( const size_t& x, const size_t& y,  const complex_t& f ) { norms[x * cols + y] = f; calculate[x*cols+y] = false; }
void Bilinears:: setZero( const size_t& x, const size_t& y ){ markAsZero[x * cols + y] = true; }
void Bilinears::resetCalculateFlags()
{ 
  for(size_t i = 0 ;i < calculate.size(); ++i) calculate[i] = true; 
} 
complex_t& Bilinears::operator()( const size_t& x, const size_t& y ){ return norms[x*cols+y];}
bool   Bilinears::isZero( const size_t& x, const size_t& y ){ return markAsZero[ x*cols+y] ; } 
bool Bilinears::workToDo( const size_t& x, const size_t& y ) const { return calculate[x*cols+y] ; } // && ! markAsZero[x*cols+y] ; }
void   Bilinears::resize( const size_t& r, const size_t& c)
{
  rows = r;
  cols = c;
  norms.resize( r * c );
  markAsZero.resize( r * c );
  calculate.resize(  r * c );
  for( size_t i = 0 ; i < r*c; ++i ){
    norms[i]      = 0;
    markAsZero[i] = false;
    calculate[i]  = true;
  }
}

void Integrator::integrateBlock()
{
  #pragma omp parallel for
  for ( size_t roll = 0; roll < N; ++roll ) {
    float_v re( 0.f );
    float_v im( 0.f );
    auto b1 = m_cache.data() + m_integrals[roll].i * m_cache.nBlocks();
    auto b2 = m_cache.data() + m_integrals[roll].j * m_cache.nBlocks();
    for ( size_t i = 0; i < m_cache.nBlocks(); ++i ) {
      auto c = b1[i] * conj(b2[i]);
      #if ENABLE_AVX2 
      re = fmadd(re, m_weight[i], real(c) );
      im = fmadd(im, m_weight[i], imag(c) );
      #else
      re = re + m_weight[i] * real(c);
      im = im + m_weight[i] * imag(c); 
      #endif
    }
    m_integrals[roll].transfer( utils::sum_elements( complex_v(re, im) ) / m_norm );
  }
  m_counter = 0;
}

bool Integrator::isReady()            const { return m_events != nullptr; }

void Integrator::queueIntegral(const size_t& c1, 
    const size_t& c2, 
    const size_t& i, 
    const size_t& j, 
    Bilinears* out, 
    const bool& sim)
{
  if( !out->workToDo(i,j) ) return;
  if( sim ) 
    addIntegralKeyed( c1, c2, [out,i,j]( const complex_t& val ){ out->set(i,j,val); if( i != j ) out->set(j,i, std::conj(val) ); } );
  else 
    addIntegralKeyed( c1, c2, [out,i,j]( const complex_t& val ){ out->set(i,j,val); } );
}

void Integrator::addIntegralKeyed( const size_t& c1, const size_t& c2, const TransferFCN& tFunc )
{
  m_integrals[m_counter++] = Integral<complex_t>(c1,c2,tFunc);
  if ( m_counter == N ) integrateBlock(); 
}

void Integrator::queueIntegral(const size_t& i, const size_t& j, complex_t* result)
{
  addIntegralKeyed(i, j, [result](const complex_t& val){ *result = val ; } ); 
}

void Integrator::flush()
{
  if ( m_counter == 0 ) return;
  integrateBlock();
}

#if ENABLE_AVX2
template <> complex_t Integrator::get( const unsigned& index, const unsigned& evt ) const 
{
  return utils::at( m_cache( evt/utils::size<float_v>::value, index), evt % utils::size<float_v>::value ); 
}
#endif

template <> complex_v Integrator::get( const unsigned& index, const unsigned& evt ) const 
{
  return m_cache(evt, index); 
}
