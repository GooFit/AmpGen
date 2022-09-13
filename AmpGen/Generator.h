#ifndef AMPGEN_GENERATOR_H
#define AMPGEN_GENERATOR_H

#include "AmpGen/EventList.h"
#if USE_SIMD 
  #include "AmpGen/EventListSIMD.h"
#endif
#include "AmpGen/simd/utils.h"
#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ProgressBar.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen
{
  template <typename phaseSpace_t = PhaseSpace>
    class Generator
    {
        #if USE_SIMD 
          using eventlist_t = EventListSIMD;
        #else
          using eventlist_t = EventList;
        #endif

      private:
        EventType    m_eventType;
        phaseSpace_t m_gps;
        size_t       m_generatorBlock = {5000000};
        TRandom*     m_rnd            = {gRandom};
        bool         m_normalise      = {true};

      public:
        template <typename... ARGS> explicit Generator( const ARGS&... args ) : m_gps(args...)
        {
          m_eventType = m_gps.eventType();
          if( m_rnd != gRandom ) setRandom( m_rnd );
        }
  
        phaseSpace_t phsp() { return m_gps; }
  
        void setRandom( TRandom* rand )
        {
          m_rnd = rand;
          m_gps.setRandom( m_rnd );
        }
        void setBlockSize( const size_t& blockSize ) { m_generatorBlock = blockSize; }
        void setNormFlag( const bool& normSetting )  { m_normalise = normSetting; }

        void fillEventListPhaseSpace( eventlist_t& events, const size_t& N)
        {
          if constexpr( std::is_same<phaseSpace_t, PhaseSpace>::value )
          {
            constexpr auto w = utils::size<float_v>::value; 
            for( unsigned i = 0 ; i != N ; ++i ){
              double* addr = reinterpret_cast<double*>( events.block( i/w  ))+ i % w;
              m_gps.fill(addr, w);     
            }     
          }
          else {
            auto it = events.begin();
            while( it != events.end() )
            {
              *it = m_gps.makeEvent();
              ++it;
            }
          }
        }
        template <typename pdf_t> double getMax(const pdf_t& pdf, const eventlist_t& events) const 
        {
          double max = 0.;
          for ( const auto& evt : events ) 
          {
            auto value             = evt.genPdf();
            if( std::isnan(value) ){ 
              ERROR( value );
              evt.print(); 
            }
            else if ( value > max ) max = value;
          }
          DEBUG( "Returning normalisation constant = " << max ); 
          return max;
        }

        template <typename eventList_t, typename pdf_t> void fillEventList( pdf_t& pdf, eventList_t& list, const size_t& N)
        {
          if ( m_rnd == nullptr ) 
          {
            ERROR( "Random generator not set!" );
            return;
          }
          double maxProb   = m_normalise ? 0 : 1;
          auto size0       = list.size();
          double totalGenerated = 0; 
          pdf.reset( true );
          ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
          ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
          std::vector<bool> efficiencyReport(m_generatorBlock, false); 

          while ( list.size() - size0 < N ) {
            eventlist_t mc( m_eventType );
            mc.resize(m_generatorBlock);
            t_phsp.start();
            fillEventListPhaseSpace(mc, m_generatorBlock);
            t_phsp.stop();
            t_eval.start();
            pdf.setEvents(mc);
            pdf.prepare();
            auto previousSize = list.size();
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for ( size_t block=0; block < mc.nBlocks(); ++block )
            {
              mc.setGenPDF(block, pdf(mc.block(block), block) / mc.genPDF(block) );
            }
            maxProb = maxProb == 0 ? 1.5 * getMax(pdf, mc) : maxProb; 
            DEBUG( "Norm: " << maxProb );            
           // if constexpr ( std::is_same<phaseSpace_t, TreePhaseSpace>::value ) m_gps.recalculate_weights(mc); 

            t_eval.stop();
            t_acceptReject.start(); 
            totalGenerated += mc.size();
            for(const auto& event : mc)
            { 
              if ( event.genPdf()  > maxProb ) {
                std::cout << std::endl; 
                WARNING( "PDF value exceeds norm value: " << event.genPdf() << " > " << maxProb );
                event.print();
              }
              if ( event.genPdf() > maxProb * m_rnd->Rndm() ){
                list.push_back(event);
                list.rbegin()->setGenPdf( pdf(event) );
                efficiencyReport[event.index()] = true; 
              }
              else efficiencyReport[event.index()] = false; 
              if ( list.size() - size0 == N ) break; 
            }
            /// if constexpr ( std::is_same<phaseSpace_t, TreePhaseSpace>::value ) maxProb = 0; 
            t_acceptReject.stop(); 
            double efficiency = 100. * ( list.size() - previousSize ) / (double)m_generatorBlock;
            pb.print( double(list.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
            if ( list.size() == previousSize ) {
              ERROR( "No events generated, PDF: " << type_string<pdf_t>() << " is likely to be malformed" );
              break;
            }
          } 
          pb.finish();
          t_total.stop();
          INFO("Generated " << N << " events in " << t_total << " ms");
          INFO("Generating phase space : " << t_phsp         << " ms"); 
          INFO("Evaluating PDF         : " << t_eval         << " ms"); 
          INFO("Accept/reject          : " << t_acceptReject << " ms"); 
          INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");
        }
        template <typename pdf_t, typename = typename std::enable_if<!std::is_integral<pdf_t>::value>::type>
          EventList generate(pdf_t& pdf, const size_t& nEvents )
          {
            eventlist_t evts( m_eventType );
            fillEventList( pdf, evts, nEvents );
            EventList output( m_eventType );
            for( const auto& event : evts ) output.emplace_back( event );
            return output; 
          }
        EventList generate(const size_t& nEvents)
        {
          eventlist_t evts( m_eventType );
          evts.resize( nEvents ); 
          fillEventListPhaseSpace( evts, nEvents);
          EventList output( m_eventType );
          for( const auto& event : evts ) output.emplace_back( event );
          return output; 
        }
    };

  template <class FCN> class PDFWrapper 
  {
    public:
      void prepare(){};
      void setEvents( AmpGen::EventList& /*evts*/ ){};
      double prob_unnormalised( const AmpGen::Event& evt ) const { return m_fcn(evt); }
      explicit PDFWrapper( const FCN& fcn ) : m_fcn(fcn) {}
      size_t size() const { return 0; }
      void reset( const bool& /*flag*/ = false ){};

    private:
      FCN m_fcn;
  };

  template <class FCN> PDFWrapper<FCN> make_pdf(const FCN& fcn){ return PDFWrapper<FCN>(fcn) ; }

  /** @function PyGenerate

    Wrapper around the a phase space generator from a stringy event type to be used with python / numpy.
    */
  extern "C" void PyGenerate( const char* eventType, double* out, const unsigned int size );
} // namespace AmpGen
#endif
