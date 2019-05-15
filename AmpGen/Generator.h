#ifndef AMPGEN_GENERATOR_H
#define AMPGEN_GENERATOR_H

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Utilities.h"

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
    Generator( const ARGS&... args )
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
    void fillEventListPhaseSpace( EventList& list, const size_t& N, const size_t& cacheSize = 0 )
    {
      fillEventListPhaseSpace( list, N, cacheSize, []( const Event& evt ) { return 1; } );
    }
    void setBlockSize( const size_t& blockSize ) { m_generatorBlock = blockSize; }
    void setNormFlag( const bool& normSetting ) { m_normalise = normSetting; }

    template <class HARD_CUT>
    void fillEventListPhaseSpace( EventList& list, const size_t& N, const size_t& cacheSize, HARD_CUT cut )
    {
      unsigned int rejected = 0;
      #ifdef DEBUGLEVEL
      auto t_start          = std::chrono::high_resolution_clock::now();
      #endif
      list.reserve( N );
      while( list.size() < N ){  
        Event newEvent = m_gps.makeEvent( cacheSize );
        newEvent.setWeight( 1 );
        if ( cut( newEvent ) ) list.push_back( newEvent );
        else rejected ++;
      }
      #ifdef DEBUGLEVEL
      auto t_end  = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration<double, std::milli>( t_end - t_start ).count();
      #endif
      DEBUG( "Stage 1 efficiency = " << 100. * list.size() / ( list.size() + rejected ) << "%, yield = " << list.size()
                                     << " time = " << time );
    }
    template <class PDF>
    void fillEventList( PDF& pdf, EventList& list, const size_t& N )
    {
      fillEventList( pdf, list, N, []( const Event& evt ) { return 1; } );
    }

    template <class PDF, class HARD_CUT>
    void fillEventList( PDF& pdf, EventList& list, const size_t& N, HARD_CUT cut )
    {
      if ( m_rnd == nullptr ) {
        ERROR( "Random generator not set!" );
        return;
      }
      double normalisationConstant = m_normalise ? 0 : 1;
      size_t size0                 = list.size();
      auto tStartTotal             = std::chrono::high_resolution_clock::now();
      pdf.reset( true );
      while ( list.size() - size0 < N ) {
        auto t_start = std::chrono::high_resolution_clock::now();
        EventList mc( m_eventType );
        fillEventListPhaseSpace( mc, m_generatorBlock, pdf.size(), cut );
 
        pdf.setEvents( mc );
        pdf.prepare();

        if ( normalisationConstant == 0 ) {
          double max = 0;
          for ( auto& evt : mc ) {
            double value           = pdf.prob_unnormalised( evt );
            if ( value > max ) max = value;
          }
          normalisationConstant = max * 1.5;
          INFO( "Setting normalisation constant = " << normalisationConstant );
        }

        unsigned int previousSize = list.size();
        for ( auto& evt : mc ) {
          double value = pdf.prob_unnormalised( evt );
          if ( value > normalisationConstant ) {
            WARNING( "PDF value exceeds norm value: " << value << " " << normalisationConstant );
          }
          if ( value > normalisationConstant * m_rnd->Rndm() ) {
            evt.setGenPdf( value );
            list.push_back( evt );
          }
          if ( list.size() - size0 == N ) break;
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        // double time = std::chrono::duration<double, std::milli>(t_end-t_stage2).count() ;
        double timeTotal = std::chrono::duration<double, std::milli>( t_end - t_start ).count();
        INFO( "Generator Efficiency = " << 100. * ( list.size() - previousSize ) / (double)m_generatorBlock
                                        << "% integrated yield = " << list.size() << ", time = " << timeTotal << "ms" );
        if ( list.size() == previousSize ) {
          ERROR( "No events generated, PDF: " << typeof<PDF>() << " is likely to be malformed" );
          break;
        }
      }
      double time =
          std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
      INFO( "Generated " << N << " events in " << time << " ms" );
    }
    template <class PDF, 
              class = typename std::enable_if<!std::is_integral<PDF>::value>::type>
    EventList generate(PDF& pdf, const size_t& nEvents )
    {
      EventList evts( m_eventType );
      fillEventList( pdf, evts, nEvents );
      return evts;
    }
   // template <class N, 
   //           class = typename std::enable_if<std::is_integral<N>::value>::type>
    EventList generate(const size_t& nEvents, const size_t& cacheSize=0)
    {
      EventList evts( m_eventType );
      fillEventListPhaseSpace( evts, nEvents, cacheSize );
      return evts;
    }
  };
  template <class FCN>
  class PDFWrapper {
    FCN m_fcn;
    public:
    void prepare(){};
    void setEvents( AmpGen::EventList& evts ){};
    double prob_unnormalised( const AmpGen::Event& evt ) const { return m_fcn(evt); }
    PDFWrapper( const FCN& fcn ) : m_fcn(fcn) {}
    size_t size() const { return 0; }
    void reset( const bool& flag = false ){};
  };
  extern "C" void PyGenerate( const char* eventType, double* out, const unsigned int size );
} // namespace AmpGen
#endif
