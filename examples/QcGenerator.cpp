#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/MinuitParameterSet.h"

#include <TFile.h>
#include <TRandom3.h>
#include <TRandom.h>
#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

using namespace AmpGen;
using namespace std::complex_literals;

class FixedLibPdf 
{
  public: 
    FixedLibPdf() = default; 
    FixedLibPdf(const EventType& type, MinuitParameterSet& mps) : 
      FixedLibPdf(NamedParameter<std::string>(type.decayDescriptor()+"::lib").getVal()) 
    {
      INFO("Constructing: " << type << " flib = " <<type.decayDescriptor()+"::lib" );
    }
    FixedLibPdf(const std::string& lib)
    {
      void* handle = dlopen( lib.c_str(), RTLD_NOW );
      if ( handle == nullptr ) ERROR( dlerror() );
      INFO("Constructing from " << lib );
      amp = AmpGen::DynamicFCN<complex_t( const double*, const int& )>( handle, "AMP" );
    }
    void prepare(){};
    void setEvents( AmpGen::EventList& evts ){};
    double prob_unnormalised( const AmpGen::Event& evt ) const { return std::norm(getValNoCache(evt)); }
    complex_t getValNoCache( const AmpGen::Event& evt ) const { return amp(evt,+1); }
    size_t size() { return 0; }
    void reset( const bool& flag = false ){};
  private:
    AmpGen::DynamicFCN<complex_t( const double*, const int& )> amp;
};

struct DTEvent 
{
  AmpGen::Event signal;
  AmpGen::Event    tag;
  double prob;
  DTEvent() : signal(0,0,0), tag(0,0,0) {};
  DTEvent( const AmpGen::Event& signal, const AmpGen::Event& tag ) : signal(signal), tag(tag) {};
  void set( const AmpGen::Event& s1, const AmpGen::Event& s2 ) { signal.set(s1); tag.set(s2); };
  void invertParity(){
    for( size_t i = 0 ; i < signal.size(); ++i ) if( i % 4 != 3 ) signal[i] *= -1;
    for( size_t i = 0 ; i < tag.size(); ++i )    if( i % 4 != 3 ) tag[i] *= -1;
  }
};

struct DTEventList : public std::vector<DTEvent> 
{
  AmpGen::EventType m_sigType; 
  AmpGen::EventType m_tagType;
  DTEventList( const AmpGen::EventType& signal, const AmpGen::EventType& tag ) : m_sigType(signal), m_tagType(tag) {}
  std::string particleName(const AmpGen::EventType& type, const size_t& j);
  TTree* tree(const std::string& name);
};

class DTYieldCalculator {
  public:
    DTYieldCalculator(const double& productionCrossSection = 3260) : 
      productionCrossSection(productionCrossSection){}

    double operator()(const double& lumi, 
        const AmpGen::EventType& t_signal, 
        const AmpGen::EventType& t_tag, 
        const bool& print = false);
    double bf( const AmpGen::EventType& type ) const;
  private: 
    double productionCrossSection;
    std::map<std::string, double> getKeyed( const std::string& name );
    std::map<std::string, double> branchingRatios = {getKeyed("BranchingRatios")};
    std::map<std::string, double> efficiencies    = {getKeyed("Efficiencies")};
};

void add_CP_conjugate( MinuitParameterSet& mps );

template <class PDF> struct normalised_pdf {
  PDF       pdf; 
  complex_t norm;
  normalised_pdf() = default; 
  normalised_pdf(const EventType& type, MinuitParameterSet& mps, const DTYieldCalculator& yc) : pdf(type, mps)
  {
    ProfileClock pc; 
    auto normEvents = Generator<PhaseSpace>(type).generate(1e6);
    double n =0;
    pdf.prepare();
    #pragma omp parallel for reduction(+:n)
    for(size_t i =0; i<normEvents.size(); ++i) 
      n += std::norm(pdf.getValNoCache(normEvents[i]));
    auto it = mps.find( type.decayDescriptor() + "::strongPhase");
    norm = sqrt(yc.bf(type)/n);
    if( it != nullptr ) norm *= exp( 1i * it->mean() * M_PI/180. );
    pc.stop();
    INFO(type << " Time to construct: " << pc << "[ms], norm = " << norm  << " " << typeof<PDF>() );
  }
  complex_t operator()(const Event& event){ return norm * pdf.getValNoCache(event); }
};

struct ModelStore {
  MinuitParameterSet*                                 mps;
  DTYieldCalculator                                   yieldCalculator;
  std::map<std::string, normalised_pdf<CoherentSum>>  genericModels;
  std::map<std::string, normalised_pdf<FixedLibPdf>>  flibModels;
  ModelStore(MinuitParameterSet* mps, const DTYieldCalculator& yc) : mps(mps), yieldCalculator(yc) {}
  template <class T> normalised_pdf<T>& get(const EventType& type, std::map<std::string, normalised_pdf<T>>& container)
  {
    auto key = type.decayDescriptor();
    if( container.count(key) == 0 ) 
      container[key] = normalised_pdf<T>(type, *mps, yieldCalculator);
    return container[key]; 
  }
  template <class T>    normalised_pdf<T>& find(const EventType& type);
};

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

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv, "Toy simulation for Quantum Correlated Ψ(3770) decays");
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  size_t hwt = std::thread::hardware_concurrency();
  size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
  double luminosity   = NamedParameter<double>("Luminosity"  , 818.3       , "Luminosity to generate. Defaults to CLEO-c integrated luminosity.");
  size_t nEvents      = NamedParameter<size_t>("nEvents"     , 0           , "Can also generate a fixed number of events per tag, if unspecified use the CLEO-c integrated luminosity.");
  size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
  bool   poissonYield = NamedParameter<bool  >("PoissonYield", true        , "Flag to include Poisson fluctuations in expected yields (only if nEvents is not specified)");
  double crossSection = NamedParameter<double>("CrossSection", 3.26 * 1000 , "Cross section for e⁺e⁻ → Ψ(3770) → DD'");
  std::string output  = NamedParameter<std::string>("Output" , "ToyMC.root", "File containing output events"); 
  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();

  gRandom = new TRandom3(seed);
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO("Setting " << nThreads << " fixed threads for OpenMP");
  omp_set_dynamic(0);  
#endif
  MinuitParameterSet MPS; 
  MPS.loadFromStream();
  add_CP_conjugate( MPS );
  EventType signalType( pNames );
  TFile* f = TFile::Open( output.c_str() ,"RECREATE");
  auto yc = DTYieldCalculator(crossSection);
  if( nEvents == 0 ) 
    INFO("Generating events using PDG/efficiencies with luminosity = " << luminosity << " pb⁻¹; σ = " << crossSection << " pb" );
  else INFO("Generating " << nEvents << " per sample");
  ModelStore models(&MPS, yc); 
  for( auto& tag : tags ){
    auto tokens       = split(tag, ' ');
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    type = tagParticle.eventType();
    double yield_noQC = yc(luminosity,signalType,type,true);
    std::string flib         = NamedParameter<std::string>( type.decayDescriptor() +"::lib", "");
    if( flib == "" ){
      auto generator    = Psi3770<CoherentSum,CoherentSum>(models, signalType, type); 
      double rho        = generator.rho();
      double yield = nEvents; 
      if( nEvents == 0 && poissonYield  ) yield = gRandom->Poisson(yield_noQC*rho);
      if( nEvents == 0 && !poissonYield ) yield = yield_noQC*rho;
      INFO( "Tag = " << type << " Expected Yield [incoherent] = " << yield_noQC << " rho = " << rho << " requested = " << yield );
      generator.generate(yield).tree(tokens[0])->Write();
    }
    else {
      auto generator    = Psi3770<CoherentSum,FixedLibPdf>(models, signalType, type); 
      double rho        = generator.rho();
      double yield = nEvents; 
      if( nEvents == 0 && poissonYield  ) yield = gRandom->Poisson(yield_noQC*rho); 
      if( nEvents == 0 && !poissonYield ) yield = yield_noQC*rho;  
      INFO( "Tag = " << type << " Expected Yield [incoherent] = " << yield_noQC << " rho = " << rho << " requested = " << yield );
      generator.generate(yield).tree(tokens[0])->Write();
    }
  }
  f->Close(); 
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );
}

void add_CP_conjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  for( auto& param : mps ){
    const std::string name = param->name();
    size_t pos=0;
    std::string new_name = name; 
    int sgn=1;
    if( name.find("::") != std::string::npos ){
      pos = name.find("::");
      auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
      if( props != 0 ) new_name = props->anti().name() + name.substr(pos); 
    }
    else { 
      auto tokens=split(name,'_');
      std::string reOrIm = *tokens.rbegin();
      std::string pname   = tokens[0];
      if ( reOrIm == "Re" || reOrIm == "Im" ){
        auto p = Particle( pname ).conj();
        sgn = reOrIm == "Re" ? p.quasiCP() : 1; 
        new_name = p.uniqueString() +"_"+reOrIm;
      }
      else if( tokens.size() == 2 ) {
        auto props = AmpGen::ParticlePropertiesList::get( name );
        if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
      }
    }
    if( mps.find( new_name ) == nullptr ){
      tmp.push_back( new MinuitParameter(new_name, Flag::Free, sgn * param->mean(), param->err(), 0, 0));
    }
  }
  for( auto& p : tmp ) mps.add( p );
}

std::map<std::string, double> DTYieldCalculator::getKeyed(const std::string& name)
{
  std::vector<std::string> things = AmpGen::NamedParameter<std::string>(name).getVector();
  std::map< std::string , double > branchingRatios; 
  for( auto& thing : things ){
    auto tokens = AmpGen::split( thing, ' ' );
    AmpGen::Particle p(tokens[0]);
    branchingRatios[p.uniqueString()]        = stod(tokens[1]);
    branchingRatios[p.conj().uniqueString()] = stod(tokens[1]);
  }
  return branchingRatios;
}

double DTYieldCalculator::operator()(const double& lumi, 
    const AmpGen::EventType& t_signal, 
    const AmpGen::EventType& t_tag, 
    const bool& print){
  auto statisticalFactor = 2;
  if( t_signal == t_tag || t_signal == t_tag.conj(false,true) ) statisticalFactor = 1;
  auto signal       = AmpGen::Particle(t_signal.decayDescriptor()).uniqueString(); 
  auto tag          = AmpGen::Particle(t_tag.decayDescriptor()).uniqueString();
  auto signalBar    = AmpGen::replaceAll( signal, "D0", "Dbar0");
  auto tagBar       = AmpGen::replaceAll( tag,    "D0", "Dbar0");
  auto eff          = [this](const std::string& tag) -> double { auto it = efficiencies.find(tag); 
    if( it == efficiencies.end() ){
      WARNING("Efficiency for final state: " << tag << " not found");
      return 1;
    }
    return it->second;
  };
  double efficiency = eff(signal) * eff(tag);
  double br         = branchingRatios[signal] * branchingRatios[tagBar] + branchingRatios[signalBar] * branchingRatios[tag];
  double totalDDbar = lumi * productionCrossSection; 
  if( print ){
    INFO("Expected yield for final state: " << t_signal << " vs " << t_tag );
    INFO("Total DDbar = " << totalDDbar );
    INFO("Efficiency  = " << efficiency );
    INFO("BR          = " << branchingRatios[signal] <<" x " << branchingRatios[tagBar] 
        << " + " << branchingRatios[signalBar] << " x " <<  branchingRatios[tag] 
        << " = " << br );
    INFO("Total       = " << statisticalFactor * totalDDbar * efficiency * br );
  }
  return statisticalFactor * totalDDbar * efficiency * br; 
}

double DTYieldCalculator::bf( const AmpGen::EventType& type ) const { 
  auto  p = AmpGen::Particle( type.decayDescriptor() );
  auto it = branchingRatios.find(p.decayDescriptor());
  if( it != branchingRatios.end() ) return it->second; 
  else {
    ERROR("Tag: " << p.decayDescriptor() << " not found in branching ratios");
    return 0;
  }
}
std::string DTEventList::particleName(const AmpGen::EventType& type, const size_t& j)
{
  auto count = type.count(j);
  if( count.second == 1 ) return programatic_name(type[j]);
  return programatic_name(type[j])+std::to_string(count.first);
}

TTree* DTEventList::tree(const std::string& name)
{
  DTEvent tmp(at(0).signal, at(0).tag);
  std::vector<int> id_sig(m_sigType.size()), 
    ids_sig(m_sigType.size()), 
    id_tag(m_tagType.size()),
    ids_tag(m_tagType.size());

  TTree* outputTree = new TTree(name.c_str(),name.c_str());
  for(size_t i = 0 ; i < m_sigType.size(); ++i )
  {
    outputTree->Branch((particleName(m_sigType, i)+"_PX").c_str(), &tmp.signal[4*i+0]); 
    outputTree->Branch((particleName(m_sigType, i)+"_PY").c_str(), &tmp.signal[4*i+1]); 
    outputTree->Branch((particleName(m_sigType, i)+"_PZ").c_str(), &tmp.signal[4*i+2]); 
    outputTree->Branch((particleName(m_sigType, i)+"_E").c_str(),  &tmp.signal[4*i+3]);
    outputTree->Branch((particleName(m_sigType, i)+"_ID").c_str(), &id_sig[i]);
    ids_sig[i] = ParticlePropertiesList::get( m_sigType[i] )->pdgID();
  }  
  for(size_t i = 0 ; i < m_tagType.size(); ++i )
  {
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PX").c_str(), &tmp.tag[4*i+0]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PY").c_str(), &tmp.tag[4*i+1]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PZ").c_str(), &tmp.tag[4*i+2]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_E").c_str(),  &tmp.tag[4*i+3]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_ID").c_str(), &id_tag[i]);
    ids_tag[i] = ParticlePropertiesList::get( m_tagType[i] )->pdgID();
  }
  for( auto& evt: *this ){
    bool swap = gRandom->Uniform() > 0.5;
    tmp.set(evt.signal, evt.tag);
    if( swap ) tmp.invertParity();
    for(size_t i=0; i != m_sigType.size(); ++i)
      id_sig[i] = swap ? -ids_sig[i] : ids_sig[i];
    for(size_t i=0; i != m_tagType.size(); ++i)
      id_tag[i] = swap ? -ids_tag[i] : ids_tag[i];
    outputTree->Fill();
  }
  return outputTree;
}    

template <> normalised_pdf<CoherentSum>& ModelStore::find(const EventType& type){ return get<CoherentSum>(type, genericModels); }
template <> normalised_pdf<FixedLibPdf>& ModelStore::find(const EventType& type){ return get<FixedLibPdf>(type, flibModels); }

