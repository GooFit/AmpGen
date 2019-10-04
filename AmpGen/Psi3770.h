
#include "AmpGen/NormalisedPdf.h"

using namespace AmpGen;
using namespace std::complex_literals;
template <class T1, class T2> class Psi3770 {
  private:
    EventType           m_signalType;
    EventType           m_tagType;
    normalised_pdf<T1>& m_signal; 
    normalised_pdf<T1>& m_signalBar; 
    normalised_pdf<T2>& m_tag; 
    normalised_pdf<T2>& m_tagBar; 
    PhaseSpace          m_signalPhsp;
    PhaseSpace          m_tagPhsp;
    PhaseSpace          m_headPhsp;
    bool                m_printed     = {false};
    bool                m_ignoreQc    = {NamedParameter<bool>("IgnoreQC",false)};
    size_t              m_blockSize   = {1000000}; 
  public:
    template <class T> std::tuple<double,double,double> 
      z(T& t1, T& t2, const EventType& type) const 
      {
        double n1(0), n2(0), zR(0), zI(0);
        auto normEvents = Generator<PhaseSpace>(type).generate(m_blockSize);
#pragma omp parallel for reduction(+:zR,zI,n1,n2)
        for(size_t i = 0; i < m_blockSize; ++i){
          auto p1 = t1(normEvents[i]);
          auto p2 = t2(normEvents[i]);
          auto f  = p2 * std::conj(p1);
          zR      += std::real(f);
          zI      += std::imag(f);
          n1      += std::norm(p1);
          n2      += std::norm(p2);
        }
        complex_t z2(zR,zI);
        auto arg = std::arg(z2);
        return std::make_tuple( std::abs(z2/sqrt(n1*n2)), 180 * ( arg > 0 ? arg : 2*M_PI+arg) /M_PI, sqrt(n1/n2) );
      } 
    Psi3770(ModelStore& models,
        const EventType& signalType, 
        const EventType& tagType) :
      m_signalType(signalType)                  ,
      m_tagType   (tagType)                     ,
      m_signal    (models.find<T1>(m_signalType)),
      m_signalBar (models.find<T1>(m_signalType.conj(true))),
      m_tag       (models.find<T2>(m_tagType))              ,
      m_tagBar    (models.find<T2>(m_tagType.conj(true)))   ,
      m_signalPhsp(m_signalType )               ,
      m_tagPhsp   (m_tagType )                  ,
      m_headPhsp  (EventType({"psi(3770)0","D0","Dbar0"}))
      { 
        INFO("Type         = " << m_signalType  << " x " << m_tagType.conj(true) << "] - [" << m_signalType.conj(true) << " x " << m_tagType << "]");
        auto z1 = z(m_signal, m_signalBar, m_signalType);
        auto z2 = z(m_tag   , m_tagBar   , m_tagType);
        INFO( "Signal: R = " << round(std::get<0>(z1),5) << "; δ = " << round(std::get<1>(z1),3) << "° ; r = " << round(std::get<2>(z1),5) );
        INFO( "Tag   : R = " << round(std::get<0>(z2),5) << "; δ = " << round(std::get<1>(z2),3) << "° ; r = " << round(std::get<2>(z2),5) );
      }
    complex_t operator()(const DTEvent& event )                 const { return operator()( event.signal, event.tag ); }
    complex_t operator()(const Event& signal, const Event& tag) const { return m_signal(signal)*m_tagBar(tag) - m_signalBar(signal) * m_tag(tag); }
    double P(const DTEvent& event)                              const { return m_ignoreQc ? prob_noQC(event) : std::norm(operator()(event)); }
    double prob_noQC (const DTEvent& event)                     const { return std::norm(m_signal(event.signal)*m_tagBar(event.tag)) + std::norm(m_signalBar(event.signal)*m_tag(event.tag)); } 
    std::complex<double> fSig(const DTEvent& event)             const {return m_signal(event.signal);}
    std::complex<double> fSigCC(const DTEvent& event)             const {return m_signalBar(event.signal);}
    std::complex<double> fTag(const DTEvent& event)             const {return m_tag(event.tag);}
    std::complex<double> fTagCC(const DTEvent& event)             const {return m_tagBar(event.tag);}
    double sigS (const DTEvent& event, std::vector<long unsigned int> ind)  const {return event.signal.s(ind);}
    double tagS (const DTEvent& event, std::vector<long unsigned int> ind)  const {return event.tag.s(ind);}

    DTEvent generatePhaseSpace()                                      { return DTEvent( m_signalPhsp.makeEvent(), m_tagPhsp.makeEvent() ); }
    DTEventList generate( const size_t& N )
    {
      DTEventList output( m_signalType, m_tagType );
      ProgressBar pb(60, trimmedString(__PRETTY_FUNCTION__));
      auto tStartTotal = std::chrono::high_resolution_clock::now();
      int currentSize  = 0;
      double norm      = -1;
      while( output.size() < N ){
        auto events = generatePHSP(m_blockSize);
        if(norm == -1 ){
          for(auto& event : events) if(event.prob > norm ) norm = event.prob;
          norm *= 1.5;
          INFO("Calculated normalisation of PDF = " << norm );
        }
        for( auto& event : events ){
          if( event.prob > norm * gRandom->Uniform() ) output.push_back( event ); 
          if( output.size() >= N ) break ; 
        }
        double efficiency = 100 * double(output.size() - currentSize )/double(m_blockSize);
        double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
        pb.print( double(output.size()) / double(N), " ε[gen] = " + mysprintf("%.2f",efficiency) + "% , " + std::to_string(int(time/1000.))  + " seconds" );
        currentSize = output.size();
      }
      auto psi_q = PhaseSpace(m_headPhsp); 
      auto beta  = [](const Event& event, const size_t&j){ return sqrt( event[4*j+0]*event[4*j+0] + event[4*j+1]*event[4*j+1] + event[4*j+2]*event[4*j+2] )/event[4*j+3] ; };
      auto p     = [](const Event& event, const size_t&j){ return std::make_tuple(event[4*j+0], event[4*j+1], event[4*j+2]); };
      for( auto& event : output ){
        auto psi_event = psi_q.makeEvent();
        boost( event.signal, p(psi_event,0), beta(psi_event,0));
        boost( event.tag   , p(psi_event,1), beta(psi_event,1));
      }
      double time = std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - tStartTotal ).count();
      pb.finish();
      INFO("Requested: " << N << " events t=" << time/1000 << "[ms]");
      return output;
    }
    DTEventList generatePHSP(const size_t& N, const bool& eval=true){
      DTEventList output( m_signalType, m_tagType );
      for(size_t x = 0 ; x < m_blockSize; ++x) output.emplace_back( generatePhaseSpace() ); 
#pragma omp parallel for
      for(size_t i = 0 ; i < m_blockSize; ++i ) output[i].prob = P(output[i]);
      return output;
    }
    double rho() {
      if( m_ignoreQc ) return 1;
      double withQC=0;
      double withoutQC =0;
      DTEventList evts = generatePHSP(m_blockSize, false); 
#pragma omp parallel for reduction(+:withQC,withoutQC)
      for( size_t x = 0; x<m_blockSize; ++x ){
        auto event   = evts[x]; 
        withQC    += P(event);
        withoutQC += prob_noQC(event);
      }
      return withQC / withoutQC ;
    }

};