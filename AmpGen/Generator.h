#ifndef AMPGEN_GENERATOR_H
#define AMPGEN_GENERATOR_H

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ProgressBar.h"
namespace AmpGen
{
  template <class PHASESPACE = PhaseSpace>
    class Generator
    {
      private:
        EventType  m_eventType;
        PHASESPACE m_gps;
        size_t     m_generatorBlock = {5000000};
        TRandom*   m_rnd            = {gRandom};
        bool       m_normalise      = {true};

      public:
        template <class... ARGS>
          explicit Generator( const ARGS&... args )
          : m_gps(args...)
          {
            m_eventType = m_gps.eventType();
            setRandom( m_rnd );
          }
        PHASESPACE phsp() { return m_gps; }
        void setRandom( TRandom* rand )
        {
          m_rnd = rand;
          m_gps.setRandom( m_rnd );
        }
        void fillEventListPhaseSpace( EventList& list, const size_t& N)
        {
          list.reserve( N );
          while( list.size() < N ){  
            Event newEvent = m_gps.makeEvent();
            newEvent.setWeight( 1 );
            newEvent.setIndex( list.size() );
            list.push_back( newEvent );
          }
        }
        void setBlockSize( const size_t& blockSize ) { m_generatorBlock = blockSize; }
        void setNormFlag( const bool& normSetting ) { m_normalise = normSetting; }

        template <class HARD_CUT> void fillEventListPhaseSpace( EventList& list, const size_t& N, HARD_CUT cut)
        {
          list.reserve( N );
          while( list.size() < N ){  
            Event newEvent = m_gps.makeEvent();
            newEvent.setWeight( 1 );
            if ( cut( newEvent ) ){
              newEvent.setIndex( list.size() );
              list.push_back( newEvent );
            }
          }
        }

        template <class PDF>
          void fillEventList( PDF& pdf, EventList& list, const size_t& N )
          {
            fillEventList( pdf, list, N, []( const Event& /*evt*/ ) { return 1; } );
          }

        template <class PDF, class HARD_CUT>
          void fillEventList( PDF& pdf, EventList& list, const size_t& N, HARD_CUT cut )
          {
            if ( m_rnd == nullptr ) {
              ERROR( "Random generator not set!" );
              return;
            }
            double maxProb   = m_normalise ? 0 : 1;
            auto size0       = list.size();
            auto tStartTotal = std::chrono::high_resolution_clock::now();
            pdf.reset( true );
            ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
            ProfileClock t_phsp, t_eval, t_acceptReject;
            std::vector<bool> efficiencyReport(m_generatorBlock,false); 
            
            while ( list.size() - size0 < N ) {
              EventList mc( m_eventType );
              t_phsp.start();
              fillEventListPhaseSpace(mc, m_generatorBlock, cut);
              t_phsp.stop();
              t_eval.start();
              pdf.setEvents( mc );
              pdf.prepare();
              t_eval.stop();
              if ( maxProb == 0 ) {
                double max = 0;
                for ( auto& evt : mc ) {
                  double value           = pdf(evt) / evt.genPdf();
                  if ( value > max ) max = value;
                }
                maxProb = max * 1.5;
                INFO( "Setting normalisation constant = " << maxProb );
              }
              auto previousSize = list.size();
              t_acceptReject.start();
              #ifdef _OPENMP
              #pragma omp parallel for
              #endif
              for ( size_t i=0; i < mc.size(); ++i ) 
                mc[i].setGenPdf(pdf(mc[i]) / mc[i].genPdf());

              for( size_t i=0; i != mc.size(); ++i )
              { 
                auto& evt = mc[i];
                if ( evt.genPdf() > maxProb ) {
                  std::cout << std::endl; 
                  WARNING( "PDF value exceeds norm value: " << evt.genPdf() << " > " << maxProb );
                }
                if ( evt.genPdf() > maxProb * m_rnd->Rndm() ){
                  list.push_back(evt);
                  efficiencyReport[i] = true; 
                }
                else efficiencyReport[i] = false; 
                if ( list.size() - size0 == N ) break;
              }
              m_gps.provideEfficiencyReport( efficiencyReport );
              t_acceptReject.stop();
              double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
              double efficiency = 100. * ( list.size() - previousSize ) / (double)m_generatorBlock;
              pb.print( double(list.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(time/1000.))  + " seconds" );
              if ( list.size() == previousSize ) {
                ERROR( "No events generated, PDF: " << typeof<PDF>() << " is likely to be malformed" );
                break;
              }
//              maxProb = 0;
            } 
            pb.finish();
            double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
            INFO( "Generated " << N << " events in " << time     << " ms" );
            INFO( "Generating phase space : " << t_phsp          << " ms"); 
            INFO( "Evaluating PDF         : " << t_eval          << " ms"); 
            INFO( "Doing accept/reject    : " << t_acceptReject  << " ms"); 
          }
        template <class PDF, class = typename std::enable_if<!std::is_integral<PDF>::value>::type>
          EventList generate(PDF& pdf, const size_t& nEvents )
          {
            EventList evts( m_eventType );
            fillEventList( pdf, evts, nEvents );
            return evts;
          }
        EventList generate(const size_t& nEvents)
        {
          EventList evts( m_eventType );
          fillEventListPhaseSpace( evts, nEvents);
          return evts;
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
