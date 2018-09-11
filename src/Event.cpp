#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"

using namespace AmpGen; 

void Event::print() const {
  for( unsigned int i = 0 ; i< m_event.size()/4. ; ++i ){
    double px = m_event[4*i+0];
    double py = m_event[4*i+1];
    double pz = m_event[4*i+2];
    double pE = m_event[4*i+3];
    double s = pE*pE - px*px -py*py -pz*pz;
    INFO( "P["<<i<<"] = {"<<px<<", " << py << ", " << pz << ", " << pE << "}, m = " << sqrt(fabs(s)) );
  }
}
double Event::pij( const size_t& i, const size_t& j) const {
  return get(i*4+3)*get(j*4+3) - 
    get(i*4+0)*get(j*4+0) -
    get(i*4+1)*get(j*4+1) -
    get(i*4+2)*get(j*4+2);
}

double Event::s( const size_t& index ) const {
  return get(4*index+3) * get(4*index+3)
    - get(4*index+0) * get(4*index+0)
    - get(4*index+1) * get(4*index+1)
    - get(4*index+2) * get(4*index+2);
}

double Event::s( const size_t& index1, const size_t& index2 ) const {
  return ( get(4*index1+3) + get(4*index2+3) )*( get(4*index1+3) + get(4*index2+3) ) - 
    ( get(4*index1+0) + get(4*index2+0) )*( get(4*index1+0) + get(4*index2+0) ) -
    ( get(4*index1+1) + get(4*index2+1) )*( get(4*index1+1) + get(4*index2+1) ) -
    ( get(4*index1+2) + get(4*index2+2) )*( get(4*index1+2) + get(4*index2+2) ) ;
}

double Event::s( const size_t& index1, const size_t& index2, const size_t& index3 ) const {
  return ( get(4*index1+3) + get(4*index2+3) + get(4*index3+3) ) *
    ( get(4*index1+3) + get(4*index2+3) + get(4*index3+3) ) - 
    ( get(4*index1+0) + get(4*index2+0) + get(4*index3+0) ) *
    ( get(4*index1+0) + get(4*index2+0) + get(4*index3+0) ) -
    ( get(4*index1+1) + get(4*index2+1) + get(4*index3+1) ) * 
    ( get(4*index1+1) + get(4*index2+1) + get(4*index3+1) ) -
    ( get(4*index1+2) + get(4*index2+2) + get(4*index3+2) ) * 
    ( get(4*index1+2) + get(4*index2+2) + get(4*index3+2) ) ;
}

double Event::s( const std::vector<size_t>& indices ) const {
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
