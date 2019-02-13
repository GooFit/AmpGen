#ifndef AMPGEN_SUMPDF_H
#define AMPGEN_SUMPDF_H

#include "AmpGen/EventList.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"

#include <tuple>
/* A SumPDF is the log-likelihood of the form
   -2*LL(event) = -2*log( Sum( i, P_i(event) ) )
   Where P_i are some probability density functions
   The sum is variadically unrolled at compile time, i.e. the wrapper
   is the same for 1..N pdfs. The unrolling should be properly inlined,
   hence N can be reasonably large with out afflicting either
   compile time or binary size.
   */

namespace AmpGen
{
  template <class... TYPES>
  class SumPDF
  {
  public:
    SumPDF( const TYPES&... _pdfs ) : m_pdfs( std::tuple<TYPES...>( _pdfs... ) ) {}
    std::tuple<TYPES...> m_pdfs;
    EventList* m_events;

    std::tuple<TYPES...> pdfs() const { return m_pdfs; }
    
    double getVal()
    {
//      ProfileClock pc;  
      double LL = 0;
      for_each( m_pdfs, []( auto& f ) { f.prepare(); } );
//      pc.stop();
//      ProfileClock eval;
      #pragma omp parallel for reduction( +: LL )
      for ( unsigned int i = 0; i < m_events->size(); ++i ) {
        auto prob = ((*this))(( *m_events)[i] );
        LL += log(prob);
      }
//      eval.stop();
//      std::cout << "t [prepare] = " << pc << "; t[eval] = " << eval << "ms" << std::endl; 
      return -2 * LL;
    }
    double operator()( const Event& evt )
    {
      double prob = 0;
      for_each( this->m_pdfs, [&prob, &evt]( auto& f ) { prob += f.prob( evt ); } );
      return prob;
    }
    void setEvents( EventList& events )
    {
      m_events = &events;
      for_each( m_pdfs, [&events]( auto& f ) { f.setEvents( events ); } );
    }
    
    template <class FCN>
    void forEach( const FCN& fcn )
    {
      for_each( m_pdfs, fcn );
    }
    std::size_t nPDFs() { return sizeof...( TYPES ); }
  };

  template <class... PDFS>
  auto make_pdf( PDFS&&... pdfs )
  {
    return SumPDF<PDFS...>( std::forward<PDFS>( pdfs )... );
  }
} // namespace AmpGen

#endif
