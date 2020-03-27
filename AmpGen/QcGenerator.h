/*
QcGenerator.h - a copy of Generator.h with correlated features instead!
The only difference with this class is the replacement of Event, EventList with 2 of each (signal and tag) therefore Generator expects a pdf
that can take 2 events as an argument - i.e. a CorrelatedSum
*/
#ifndef AMPGEN_QCGENERATOR_H
#define AMPGEN_QCGENERATOR_H

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ProgressBar.h"
#include "AmpGen/corrEventList.h"
namespace AmpGen
{
  template <class PHASESPACE = PhaseSpace>
    class QcGenerator
    {
      private:
        EventType  m_sigType;
        EventType  m_tagType;
        PHASESPACE m_gpsSig;
        PHASESPACE m_gpsTag;

        size_t     m_generatorBlock = {5000000};
        TRandom*   m_rnd            = {gRandom};
        bool       m_normalise      = {true};

      public:
        template <class... ARGS>
          explicit QcGenerator( PHASESPACE sigPH , PHASESPACE tagPH)
          {
            m_gpsSig = sigPH;
            m_gpsTag = tagPH;
            m_sigType = m_gpsSig.eventType();
            m_tagType = m_gpsTag.eventType();
            setRandom( m_rnd );
          }
        PHASESPACE phspSig() { return m_gpsSig; }
        PHASESPACE phspTag() { return m_gpsTag; }
        void setRandom( TRandom* rand )
        {
          m_rnd = rand;
          m_gpsSig.setRandom( m_rnd );
          m_gpsTag.setRandom( m_rnd );
        }
        void fillEventListPhaseSpace( EventList& listSig, EventList& listTag, const size_t& N, const size_t& cacheSize = 0 )
        {
          fillEventListPhaseSpace( listSig, listTag, N, cacheSize, []( const Event& evt1, const Event& evt2 ) { return 1; } );
        }
        void setBlockSize( const size_t& blockSize ) { m_generatorBlock = blockSize; }
        void setNormFlag( const bool& normSetting ) { m_normalise = normSetting; }

        template <class HARD_CUT> void fillEventListPhaseSpace( EventList& listSig, EventList& listTag, const size_t& N, const size_t& cacheSize, HARD_CUT cut )
        {
          listSig.reserve( N );
          listTag.reserve( N );
          while( listSig.size() < N ){  
            Event newEventSig = m_gpsSig.makeEvent( cacheSize );
            Event newEventTag = m_gpsTag.makeEvent( cacheSize );
            newEventSig.setWeight( 1 );
            newEventTag.setWeight( 1 );
            if ( cut( newEventSig, newEventTag ) ) listSig.push_back( newEventSig );
            if ( cut( newEventSig, newEventTag ) ) listTag.push_back( newEventTag );
          }
        }

        /// Fill EventList with a "default cut" - which lets anything through - 
        template <class PDF>
          void fillEventList( PDF& pdf, EventList& listSig, EventList& listTag, const size_t& N )
          {
            fillEventList( pdf, listSig, listTag, N, []( const Event& evt1, const Event& evt2 ) { return 1; } );
          }

        template <class PDF, class HARD_CUT>
          void fillEventList( PDF& pdf, EventList& listSig, EventList& listTag, const size_t& N, HARD_CUT cut )
          {
            if ( m_rnd == nullptr ) {
              ERROR( "Random generator not set!" );
              return;
            }
            INFO("Starting to Generate events");
            double normalisationConstant = m_normalise ? 0 : 1;
            auto size0                 = listSig.size();
            auto tStartTotal             = std::chrono::high_resolution_clock::now();
            pdf.reset( true );
            ProgressBar pb(60, trimmedString(__PRETTY_FUNCTION__) );
            ProfileClock t_phsp, t_eval, t_acceptReject;
      //      INFO("size0 = "<<size0);
    //        INFO("N = "<<N);
  //          INFO("listSig.size() = "<<listSig.size());
            while ( listSig.size() - size0 < N ) {
              EventList mcSig( m_sigType );
              EventList mcTag( m_tagType );
              EventList mcSig100( m_sigType );
              EventList mcTag100( m_tagType );
              t_phsp.start();
//              INFO("Filling Phase Space");
              fillEventListPhaseSpace( mcSig, mcTag, m_generatorBlock, pdf.size(), cut );
              //fillEventListPhaseSpace( mcSig100, mcTag100, 100 * m_generatorBlock, pdf.size(), cut );
              t_phsp.stop();
              t_eval.start();
              //INFO("Preparing CorrelatedSum");
              pdf.setEvents( mcSig, mcTag );
              pdf.setMC( mcSig, mcTag );
              pdf.prepare();

//              INFO("Prepared CorrelatedSum");
              t_eval.stop();
              if ( normalisationConstant == 0 ) {
                double max = 0;
                for ( size_t i=0; i < mcSig.size(); i++ ) {
  //                INFO("Getting value of CS");
                  double value           = pdf.prob_unnormalised( mcSig[i], mcTag[i] );
    //              INFO("Value = "<<value);
                  if ( value > max ) max = value;
                }
                normalisationConstant = max * 15;
                INFO( "Setting normalisation constant = " << normalisationConstant );
              }
              auto previousSize = listSig.size();
              t_acceptReject.start();
#ifdef _OPENMP
#pragma omp parallel for
#endif
              for ( size_t i=0;i< mcSig.size(); ++i ) mcSig[i].setGenPdf(pdf.prob_unnormalised(mcSig[i], mcTag[i]));
              for ( size_t i=0;i< mcTag.size(); ++i ) mcTag[i].setGenPdf(pdf.prob_unnormalised(mcSig[i], mcTag[i]));

              for( size_t i=0; i < mcSig.size(); ++i ){
                if ( mcSig[i].genPdf() > normalisationConstant  || mcTag[i].genPdf() > normalisationConstant ) {
                  std::cout << std::endl; 
                  if (mcSig[i].genPdf() > normalisationConstant ) WARNING( "PDF value exceeds norm value: " << mcSig[i].genPdf() << " > " << normalisationConstant );
                  if (mcTag[i].genPdf() > normalisationConstant ) WARNING( "PDF value exceeds norm value: " << mcTag[i].genPdf() << " > " << normalisationConstant );
                }
                auto threshold = normalisationConstant * m_rnd->Rndm();
                if ( mcSig[i].genPdf() > threshold ) listSig.push_back( mcSig[i] );
                if ( mcTag[i].genPdf() > threshold ) listTag.push_back( mcTag[i] );
                if ( listSig.size() - size0 == N  && listTag.size() - size0 == N) break;
              }
              t_acceptReject.stop();
              double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
              double efficiency = 100. * ( listSig.size() - previousSize ) / (double)m_generatorBlock;
              pb.print( double(listSig.size()) / double(N), " ε[gen] = " + mysprintf("%.2f",efficiency) + "% , " + std::to_string(int(time/1000.))  + " seconds" );
              if ( listSig.size() == previousSize ) {
                ERROR( "No events generated, PDF: " << typeof<PDF>() << " is likely to be malformed" );
                break;
              }
            } 
            pb.finish();
            double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
            INFO( "Generated " << N << " events in " << time     << " ms" );
            INFO( "Generating phase space : " << t_phsp          << " ms"); 
            INFO( "Evaluating PDF         : " << t_eval          << " ms"); 
            INFO( "Doing accept/reject    : " << t_acceptReject  << " ms"); 
          }
template <class PDF>
          void filterEventList( PDF& pdf, EventList& listSig, EventList& listTag, EventList& inSig, EventList& inTag, const size_t& N )
          {
            if ( m_rnd == nullptr ) {
              ERROR( "Random generator not set!" );
              return;
            }
            INFO("Starting to Generate events");
            double normalisationConstant = m_normalise ? 0 : 1;
            auto size0                 = listSig.size();
            auto tStartTotal             = std::chrono::high_resolution_clock::now();
            pdf.reset( true );
            ProgressBar pb(60, trimmedString(__PRETTY_FUNCTION__) );
            ProfileClock t_phsp, t_eval, t_acceptReject;
            while ( listSig.size() - size0 < N ) {
              t_eval.start();
             
              pdf.setEvents( inSig, inTag );
              pdf.setMC( inSig, inTag );
              pdf.prepare();
            //Normalisation constant = max(pdf)  * Norm(pdf) * 15
              t_eval.stop();
              if ( normalisationConstant == 0 ) {
                double max = 0;
                for ( size_t i=0; i < inSig.size(); i++ ) {
//                  double value           = pdf.prob_unnormalised( inSig[i], inTag[i] );
 //                 auto vals = pdf.getVals(inSig[i], inTag[i]);
  //                auto value = std::norm(vals[0]) * std::norm(vals[1]);
                  double value           = pdf.prob_unnormalised( inSig[i], inTag[i] );
                  if ( value > max ) max = value;
                }
                normalisationConstant = max * 15;
                INFO( "Setting normalisation constant = " << normalisationConstant );
              }
              auto previousSize = listSig.size();
              t_acceptReject.start();
#ifdef _OPENMP
#pragma omp parallel for
#endif
              for ( size_t i=0;i< inSig.size(); ++i ) inSig[i].setGenPdf(pdf.prob_unnormalised(inSig[i], inTag[i]));
              for ( size_t i=0;i< inTag.size(); ++i ) inTag[i].setGenPdf(pdf.prob_unnormalised(inSig[i], inTag[i]));

              for( size_t i=0; i < inSig.size(); ++i ){
              
                if ( inSig[i].genPdf() > normalisationConstant  || inTag[i].genPdf() > normalisationConstant ) {
                  std::cout << std::endl; 
                  if (inSig[i].genPdf() > normalisationConstant ) WARNING( "PDF value exceeds norm value: " << inSig[i].genPdf() << " > " << normalisationConstant );
                  if (inTag[i].genPdf() > normalisationConstant ) WARNING( "PDF value exceeds norm value: " << inTag[i].genPdf() << " > " << normalisationConstant );
                }
                auto threshold = normalisationConstant * m_rnd->Rndm();
                if ( inSig[i].genPdf() > threshold &&  inTag[i].genPdf() > threshold   )
                { listSig.push_back( inSig[i] );
               listTag.push_back( inTag[i] );
                }
              
                if ( listSig.size() - size0 == N  && listTag.size() - size0 == N) break;
              }
              t_acceptReject.stop();
              double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
              double efficiency = 100. * ( listSig.size() - previousSize ) / (double)inSig.size();
              pb.print( double(listSig.size()) / double(N), " ε[gen] = " + mysprintf("%.2f",efficiency) + "% , " + std::to_string(int(time/1000.))  + " seconds" );
              if ( listSig.size() == previousSize ) {
                ERROR( "No events generated, PDF: " << typeof<PDF>() << " is likely to be malformed" );
                break;
              }
             
            }
            pb.finish();
            double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
            INFO( "Generated " << N << " events in " << time     << " ms" );
            INFO( "Generating phase space : " << t_phsp          << " ms"); 
            INFO( "Evaluating PDF         : " << t_eval          << " ms"); 
            INFO( "Doing accept/reject    : " << t_acceptReject  << " ms"); 
          }
 
        template <class PDF, class = typename std::enable_if<!std::is_integral<PDF>::value>::type>
          corrEventList generate(PDF& pdf, const size_t& nEvents )
          {
            EventList evtsSig( m_sigType );
            EventList evtsTag( m_tagType );
            INFO("Generating event list");
            fillEventList( pdf, evtsSig, evtsTag, nEvents );
            auto cL = corrEventList(evtsSig, evtsTag);
            return cL;
          }
        corrEventList generate(const size_t& nEvents, const size_t& cacheSize=0)
        {
          EventList evtsSig( m_sigType );
          EventList evtsTag( m_tagType );
          fillEventListPhaseSpace( evtsSig, evtsTag, nEvents, cacheSize );
          auto cL = corrEventList(evtsSig, evtsTag);
          return cL;
        }

      


    };

  template <class FCN> class corrPDFWrapper 
  {
    public:
      void prepare(){};
      void setEvents( AmpGen::EventList& , AmpGen::EventList&/*evts*/ ){};
      double prob_unnormalised( const AmpGen::Event& evt1, const AmpGen::Event& evt2 ) const { return m_fcn(evt1, evt2); }
      explicit corrPDFWrapper( const FCN& fcn ) : m_fcn(fcn) {}
      size_t size() const { return 0; }
      void reset( const bool& /*flag*/ = false ){};

    private:
      FCN m_fcn;
  };

  template <class FCN> corrPDFWrapper<FCN> make_corrpdf(const FCN& fcn){ return corrPDFWrapper<FCN>(fcn) ; }

  /** @function PyGenerate

    Wrapper around the a phase space generator from a stringy event type to be used with python / numpy.
   */
  extern "C" void PyGenerate( const char* eventType, double* out, const unsigned int size );
} // namespace AmpGen
#endif
