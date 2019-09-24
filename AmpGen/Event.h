#ifndef AMPGEN_EVENT_H
#define AMPGEN_EVENT_H

#include <vector>
#include <complex>
#include <cstring>
#include <cstddef>
#include <array>
#include "AmpGen/Types.h"

namespace AmpGen { 
  
  /** @class Event 
      @brief Encapsulates the final state particles of a single event 
      
      Encapsulates the final state particles of a single event, or candidate in the language of proton-proton collisions. Typically will store (i) the event kinematics, i.e. four-momenta, (ii). a cache of complex numbers that contain intermediate calculations of the amplitude, (iii). the weight of the given event/candidate, (iv). The probability that the event was generated with, in the case of a simulated event */

  class Event {
    public:

      Event( const unsigned& N, const unsigned& cacheSize=0 );
      Event( const real_t* data, const unsigned& N, const unsigned& cacheSize=0);

      void set( const unsigned& i, const std::vector<real_t>& p );
      void set( const unsigned& i, const real_t* p );
      void set( const real_t* evt );
      void set( const unsigned& i, const real_t& p ) ;
      void swap( const unsigned int& i , const unsigned int& j );
      void setCache(const complex_t& value, const unsigned& pos) ;
      template <unsigned N> void setCache( const std::array<complex_t,N>& value, const unsigned& pos )
      {
        std::memmove( m_cache.data() + pos, value.data(), sizeof(std::array<complex_t,N>) );
      }
      void setCache( const std::vector<complex_t>& value, const unsigned& pos );
      void resizeCache( const unsigned int& new_size );
      
      unsigned   size()                                   const { return m_event.size(); } 

      real_t* pWeight()                                         { return &(m_weight); }
      real_t* pGenPdf()                                         { return &m_genPdf; }
      const real_t* address(const unsigned& ref=0)        const { return &(m_event[ref]); }
      real_t*       address(const unsigned& ref=0)              { return &(m_event[ref]); }

      unsigned cacheSize()                                const { return m_cache.size(); } 
      real_t weight()                                     const { return m_weight; } 
      real_t genPdf()                                     const { return m_genPdf; }
      real_t  operator[](const unsigned& i)               const { return m_event[i]; }
      real_t& operator[](const unsigned& i)                     { return m_event[i]; }
      operator const real_t*()                            const { return &(m_event[0]); }
      operator       real_t*()                                  { return &(m_event[0]); }

      const complex_t& getCache(const unsigned& pos)      const { return m_cache[pos]; }
      const complex_t* getCachePtr(const unsigned& pos=0) const { return &(m_cache[0]) + pos; }

      void setWeight( const real_t& weight ){ m_weight = weight ; } 
      void setGenPdf( const real_t& genPdf ){ m_genPdf = genPdf ; } 
      void extendEvent(const real_t& value) { m_event.push_back( value ); } 

      void print()      const;
      void printCache() const;

      real_t s( const unsigned& index) const ;  
      real_t s( const unsigned& index1, const unsigned& index2 ) const ;
      real_t s( const unsigned& index1, const unsigned& index2, const unsigned& index3 ) const;
      real_t s( const std::vector<unsigned>& indices ) const ;
    private:
      std::vector<real_t>    m_event; 
      std::vector<complex_t> m_cache;
      real_t                 m_genPdf = {1};
      real_t                 m_weight = {1}; 

      inline real_t get(const unsigned& index ) const { return m_event[index]; };
  };
}

#endif 
