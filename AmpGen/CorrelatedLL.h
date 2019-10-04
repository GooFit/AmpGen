#ifndef CORRELATED_LL
#define CORRELATED_LL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"

#include <tuple>

namespace AmpGen
{
  class EventList; 
  /** @class SumPDF
      @brief A pdf that contains one or more terms.

      A pdf with a probability of the form 
      @f[
        P(\psi, \phi) = \sum_{i} \sum_{j} \mathcal{P}_i (\psi) \mathcal{P}_j (\phi),
      @f] 
      where @f$ \mathcal{P}_j(\psi) @f$ are some normalised probability density functions 
      as a function of position in the phase space @f$ \psi @f$ 
      , and the sum is over the different terms, typically a signal term and then a number of background terms. 
      The pdf is also equipped with a log-likelihood of the form:
      @f[ 
        -2 \mathcal{L} = - 2 \sum_{i,j} \log \left( P_{ij}\left(\psi,\phi\right) \right)
      @f]
      and the sum over @f$ i @f$ is over some dataset of signal @f$i @f$ and tag @f$ j @$f.
      This combined functionality is largely historical and the two roles should be separated at some point in the future. 
      The sum is variadically unrolled at compile time, i.e. the wrapper
      is the same for 1..N pdfs. The unrolling should be properly inlined,
      hence N can be reasonably large with out afflicting either
      compile time or binary size. It isn't primarily used as PDF, as its primary function is 
      as a likelihood via function getVal(). 
      Typically constructed using either the make_pdf helper function or make_likelihood helper function.  */
  template <class eventListType, class... pdfTypes>
  class CorrelatedLL
  {
  private:
    typedef typename eventListType::value_type eventValueType; ///< The value type stored in the eventListType
    std::tuple<pdfTypes...> m_pdfs;                            ///< The tuple of probability density functions
    eventListType*          m_events1;                          ///< The event list to evaluate likelihoods on
    eventListType*          m_events2;                          ///< The event list to evaluate likelihoods on
    bool                    m_debug = false;

  public:
    /// Default Constructor
    CorrelatedLL() = default; 

    /// Constructor from a set of PDF functions
    CorrelatedLL( const pdfTypes&... pdfs ) : 
      m_pdfs( std::tuple<pdfTypes...>( pdfs... ) ),
      m_debug(NamedParameter<bool>("CorrelatedLL::Debug", false, "Print debug messages for CorrelatedLL")) {}

    /// Returns negative twice the log-likelihood for this PDF and the given dataset.     
    double getVal(){
        double LL = 0;
        if (m_debug) INFO("Begin getVal()");
        for_each( m_pdfs, []( auto& f ) { 
          f.prepare(); 
          //f.debugNorm();
          });

        if (m_debug) INFO("Events 1  = "<<&m_events1); 
        if (m_debug) INFO("Events 1 size  = "<<m_events1->size());
        if (m_debug) INFO("Events 2  = "<<&m_events2); 
        if (m_debug) INFO("Events 2 size  = "<<m_events2->size());
        #pragma omp parallel for reduction( +: LL )
        for ( unsigned int i = 0; i < m_events1->size(); ++i ) {
            auto prob = ((*this))( (*m_events1)[i], (*m_events2)[i]);
            //if (m_debug) INFO("Prob (from LL) = "<<prob);
            LL += log(prob);
        }
        LL *= -2;
        if (m_debug) INFO("LL = "<<LL);
        return LL;
    }

    double operator() (const Event& event1, const Event& event2){
        double prob=0;
        double prob_norm=0;
        for_each( this->m_pdfs, [&prob, &event1, &event2]( auto& f ) { prob += f.prob( event1, event2 ); } );
        for_each( this->m_pdfs, [&prob_norm, &event1, &event2]( auto& f ) { prob_norm += f.prob( event1, event2 )/f.norm(); } );
        //if (m_debug) INFO("prob_norm = "<<prob_norm);
        return prob;
    }

    void setEvents(eventListType& events1, eventListType& events2){
        m_events1 = &events1;
        m_events2 = &events2;
        for_each( m_pdfs, [&events1, &events2]( auto& f ) { f.setEvents( events1, events2 ); } );

    }
    std::tuple<pdfTypes...> pdfs() const { return m_pdfs; }

  };

    /*
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
  #endif
    */
  /** @function make_LL 

    Usage is 
    \code{cpp}
      auto pdf = make_pdf( signal, bkg1, ... );
    \endcode
    which returns a PDF that is the sum of the signal and bkg1 etc. The sum is also equipped with a likelihood, which can be used by setting 
    the data for the PDF.
    Therefore, named CorrelatedLL, it is useful to use this function to get a likelihood for a PDF containing a single term (i.e. signal or background only).  
    */
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_LL( pdfTypes&&... pdfs )
  {
    //return CorrelatedLL<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    return CorrelatedLL<eventListType, pdfTypes...>( pdfs... );
  }
  
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_likelihood( eventListType& events1, eventListType& events2, pdfTypes&&... pdfs )
  {
    auto rt = CorrelatedLL<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    rt.setEvents(events1, events2);
    return rt; 
  }
} // namespace AmpGen
#endif