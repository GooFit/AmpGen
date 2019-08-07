#include <math.h>
#include <stddef.h>
#include <complex>
#include <ostream>
#include <vector>

#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Types.h"

using namespace AmpGen; 

Event::Event( const size_t& N, const size_t& cacheSize) : 
  m_event(N), 
  m_cache(cacheSize) {
  }

Event::Event( const real_t* data, const size_t& N, const size_t& cacheSize) :
  m_event(data, data+N),
  m_cache(cacheSize) {
  }

void Event::print() const {
  size_t nParticles = m_event.size()/4;
  for( size_t i = 0 ; i < nParticles; ++i ){
    double px = m_event[4*i+0];
    double py = m_event[4*i+1];
    double pz = m_event[4*i+2];
    double pE = m_event[4*i+3];
    double s = pE*pE - px*px -py*py -pz*pz;
    INFO( "P["<<i<<"] = {"<<px<<", " << py << ", " << pz << ", " << pE << "}, m = " << sqrt(fabs(s)) );
  }
  for( size_t i = 4*nParticles; i != m_event.size(); ++i ){
    INFO( m_event[i] );
  }
}
real_t Event::s( const size_t& index ) const {
  return get(4*index+3) * get(4*index+3)
    - get(4*index+0) * get(4*index+0)
    - get(4*index+1) * get(4*index+1)
    - get(4*index+2) * get(4*index+2);
}

real_t Event::s( const size_t& index1, const size_t& index2 ) const {
  return ( get(4*index1+3) + get(4*index2+3) )*( get(4*index1+3) + get(4*index2+3) ) - 
    ( get(4*index1+0) + get(4*index2+0) )*( get(4*index1+0) + get(4*index2+0) ) -
    ( get(4*index1+1) + get(4*index2+1) )*( get(4*index1+1) + get(4*index2+1) ) -
    ( get(4*index1+2) + get(4*index2+2) )*( get(4*index1+2) + get(4*index2+2) ) ;
}

real_t Event::s( const size_t& index1, const size_t& index2, const size_t& index3 ) const {
  return ( get(4*index1+3) + get(4*index2+3) + get(4*index3+3) ) *
    ( get(4*index1+3) + get(4*index2+3) + get(4*index3+3) ) - 
    ( get(4*index1+0) + get(4*index2+0) + get(4*index3+0) ) *
    ( get(4*index1+0) + get(4*index2+0) + get(4*index3+0) ) -
    ( get(4*index1+1) + get(4*index2+1) + get(4*index3+1) ) * 
    ( get(4*index1+1) + get(4*index2+1) + get(4*index3+1) ) -
    ( get(4*index1+2) + get(4*index2+2) + get(4*index3+2) ) * 
    ( get(4*index1+2) + get(4*index2+2) + get(4*index3+2) ) ;
}

real_t Event::s( const std::vector<size_t>& indices ) const {
  if( indices.size() == 2 ) return s( indices[0], indices[1] );
  if( indices.size() == 3 ) return s( indices[0], indices[1], indices[2] );
  double E=0;        double px=0;
  double py=0;
  double pz=0;
  for( auto& i : indices ){
    E +=get(i*4+3);
    px+=get(i*4+0);
    py+=get(i*4+1);
    pz+=get(i*4+2);
  }
  return E*E -px*px - py*py - pz*pz;
}
void Event::printCache() const {
  for( unsigned int i = 0 ; i < m_cache.size(); ++i){
    INFO("Cache adddress [" << i << "] = " << m_cache[i] );
  }
}

void Event::set( const size_t& i, const std::vector<real_t>& p ){
  for( size_t j = 0 ; j < 4; ++j) m_event[4*i + j ] = p[j];
}
void Event::set( const size_t& i, const real_t* p ){
  for( size_t j = 0 ; j < 4; ++j) m_event[4*i + j ] = p[j];
}
void Event::set( const real_t* evt ){
  for( size_t i = 0 ; i < m_event.size(); ++i ) m_event[i] = *(evt + i );
}
void Event::set( const size_t& i, const real_t& p ){ m_event[i] = p ; } 

void Event::swap( const unsigned int& i , const unsigned int& j )
{
  double tmp[4];
  std::memmove( tmp, &m_event[4*j], sizeof(tmp)); 
  std::memmove( &m_event[4*j], &m_event[4*i],sizeof(tmp));
  std::memmove( &m_event[4*i], &tmp,sizeof(tmp));
}

void Event::setCache(const complex_t& value, const size_t& pos){ m_cache[pos] = value; }
void Event::setCache( const std::vector<complex_t>& value, const size_t& pos )
{
  std::memmove( m_cache.data() + pos, value.data(), sizeof(complex_t) * value.size() );
}

void Event::resizeCache( const unsigned int& new_size ){ m_cache.resize(new_size); }
