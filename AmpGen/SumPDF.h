#ifndef AMPGEN_SUMPDF_H
#define AMPGEN_SUMPDF_H

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"

#include <tuple>

namespace AmpGen
{
  class EventList; 
  /** @class SumPDF
      @brief A pdf that contains one or more terms.

      A pdf with a probability of the form 
      @f[
        P(\psi) = \sum_{j} \mathcal{P}_j (\psi),
      @f] 
      where @f$ \mathcal{P}_j(\psi) @f$ are some normalised probability density functions 
      as a function of position in the phase space @f$ \psi @f$ 
      , and the sum is over the different terms, typically a signal term and then a number of background terms. 
      The pdf is also equipped with a log-likelihood of the form:
      @f[ 
        -2 \mathcal{L} = - 2 \sum_{i} \log \left( \sum_{j} \mathcal{P}_j \left(\psi_i\right) \right)
      @f]
      and the sum over @f$ i @f$ is over some dataset. 
      This combined functionality is largely historical and the two roles should be separated at some point in the future. 
      The sum is variadically unrolled at compile time, i.e. the wrapper
      is the same for 1..N pdfs. The unrolling should be properly inlined,
      hence N can be reasonably large with out afflicting either
      compile time or binary size. It isn't primarily used as PDF, as its primary function is 
      as a likelihood via function getVal(). 
      Typically constructed using either the make_pdf helper function or make_likelihood helper function.  */
  template <class eventListType, class... pdfTypes>
  class SumPDF
  {
  private:
    typedef typename eventListType::value_type eventValueType; ///< The value type stored in the eventListType
    std::tuple<pdfTypes...> m_pdfs;                            ///< The tuple of probability density functions
    eventListType*          m_events;                          ///< The event list to evaluate likelihoods on

  public:
    /// Default Constructor
    SumPDF() = default; 

    /// Constructor from a set of PDF functions
    SumPDF( const pdfTypes&... pdfs ) : m_pdfs( std::tuple<pdfTypes...>( pdfs... ) ) {}

    /// Returns negative twice the log-likelihood for this PDF and the given dataset.     
    double getVal()
    {
      double LL = 0;
      for_each( m_pdfs, []( auto& f ) { f.prepare(); } );
      #pragma omp parallel for reduction( +: LL )
      for ( unsigned int i = 0; i < m_events->size(); ++i ) {
        auto prob = ((*this))(( *m_events)[i] );
        LL += log(prob);
      }
      return -2 * LL;
    }
    
    /// Returns the probability for the given event. 
    double operator()( const eventValueType& evt )
    {
      double prob = 0;
      for_each( this->m_pdfs, [&prob, &evt]( auto& f ) { prob += f(evt); } );
      return prob;
    }

    /// Sets the events to be summed over in the likelihood
    void setEvents( eventListType& events )
    {
      m_events = &events;
      for_each( m_pdfs, [&events]( auto& f ) { f.setEvents( events ); } );
    }
    
    /// Returns the number of PDFs contained by this function 
    std::size_t nPDFs() const { return sizeof...(pdfTypes); }

    /// Returns the tuple of PDFs used by this function
    std::tuple<pdfTypes...> pdfs() const { return m_pdfs; }
  };

  /** @function make_pdf 

    Usage is 
    \code{cpp}
      auto pdf = make_pdf( signal, bkg1, ... );
    \endcode
    which returns a PDF that is the sum of the signal and bkg1 etc. The sum is also equipped with a likelihood, which can be used by setting 
    the data for the PDF.
    Therefore, named SumPDF, it is useful to use this function to get a likelihood for a PDF containing a single term (i.e. signal or background only).  
    */
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_pdf( pdfTypes&&... pdfs )
  {
    //return SumPDF<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    return SumPDF<eventListType, pdfTypes...>( pdfs... );
  }
  
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_likelihood( eventListType& events, pdfTypes&&... pdfs )
  {
    auto rt = SumPDF<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    rt.setEvents(events);
    return rt; 
  }
} // namespace AmpGen

#endif
