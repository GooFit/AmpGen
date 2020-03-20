#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>

//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

std::vector<std::string> makeBranches(EventType Type, std::string prefix);
//template <typename PDF>
//FitResult* Fit( PDF&& pdf, pCorrelatedSum& cs, std::map<std::string, EventList> events, MinuitParameterSet& MPS , std::string logFile);
template <typename PDF>
FitResult* Fit( PDF&& ll, pCorrelatedSum& cs, EventList data1, EventList data2, EventList mc1, EventList mc2, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits );


MinuitParameterSet *  copyMPS(MinuitParameterSet& mps);

std::map<std::string, double>  getPulls(std::map<std::string, std::vector<double> > fit, std::map<std::string, std::vector<double> > init);

EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile, std::string intFile);
void writePulls(std::string fileName, std::map<std::string, double>);
std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps);


std::map<std::string, EventType> makeEventTypes(std::vector<std::string> sigName, std::string tagName);

pCorrelatedSum makeCs(std::map<std::string, EventList> events, MinuitParameterSet& MPS);
CorrelatedLL<EventList, pCorrelatedSum&> tagLL(std::map<std::string, EventList> events, MinuitParameterSet& MPS);
void add_CP_conjugate( MinuitParameterSet& mps );


template <typename pdftype>
void doPlots(std::string plotFile,  pdftype&& cs,  std::map<std::string, EventList> events, int nBins);
/*template <class Container>
void split5(const std::string& str, Container& cont,
              const std::string& delims = " ")
{
    boost::split(cont, str, boost::is_any_of(delims));
}
*/
struct Moment {
  double x;
  double xx;
  double N;
  std::vector<double> values;
  Moment() : x( 0 ), xx( 0 ), N( 0 ) {}
  void add( const double& value )
  {
    x += value;
    xx += value * value;
    N++;
    values.push_back( value );
  }
  void rescale( const double& val )
  {
    x *= val;
    xx *= ( val * val );
  }
  double val() { return x; }
  double var() { return N == 0 ? 0 : xx; }
};

int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );
  //OptionsParser::setArgs( argc, argv, "Toy simulation for Quantum Correlated Ψ(3770) decays");
  /* */
  //auto time_wall = std::chrono::high_resolution_clock::now();
  //auto time      = std::clock();
  size_t hwt = std::thread::hardware_concurrency();
  size_t nThreads     = NamedParameter<size_t>("nCores"      , hwt         , "Number of threads to use");
  //double luminosity   = NamedParameter<double>("Luminosity"  , 818.3       , "Luminosity to generate. Defaults to CLEO-c integrated luminosity.");
  //size_t nEvents      = NamedParameter<size_t>("nEvents"     , 0           , "Can also generate a fixed number of events per tag, if unspecified use the CLEO-c integrated luminosity.");
  size_t seed         = NamedParameter<size_t>("Seed"        , 0           , "Random seed to use.");
  //bool   poissonYield = NamedParameter<bool  >("PoissonYield", true        , "Flag to include Poisson fluctuations in expected yields (only if nEvents is not specified)");
  //bool   noQCFit      = NamedParameter<bool  >("noQCFit"     , false       , "Treat Signal and Tag as uncorrelated and fit the data individually");
  double crossSection = NamedParameter<double>("CrossSection", 3.26 * 1000 , "Cross section for e⁺e⁻ → Ψ(3770) → DD'");
  std::string output  = NamedParameter<std::string>("Output" , "ToyMC.root", "File containing output events"); 
  auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  auto tags           = NamedParameter<std::string>("TagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
  bool m_debug        = NamedParameter<bool>("Debug", false, "Debug QcFitter output");
    bool doDebugNorm  = NamedParameter<bool>("doDebugNorm", false, "Debug the normalisation of the pdf");
    int nBins = NamedParameter<int>("nBins", 100, "number of bins for projection");
    int nFits = NamedParameter<int>("nFits", 4, "number of repeats of mini.doFits() for debug purposes!");
    bool doProjections = NamedParameter<bool>("doProjections", true);
    bool doPCorrSum = NamedParameter<bool>("doPCorrSum", false);
//    bool doCombFit = NamedParameter<bool>("doCombFit", false, "Do a combined fit of 3 tags - at the moment this is hard coded for now");


  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
  bool doCombFit      = NamedParameter<bool>("doCombFit", false, "Do combined fit");
  bool QcGen2 = NamedParameter<bool>("QcGen2", false, "internal boolean - for new QcGenerator");
  bool doFit = NamedParameter<bool>("doFit", true, "Do the fit");
  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);
  if (intFile == ""){

  }
  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  if (m_debug) INFO("LogFile: " << logFile << "; Plots: " << plotFile );
   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  if (m_debug) INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  std::vector<std::string> varNames = {"E", "PX", "PY", "PZ"};
  //auto yc = DTYieldCalculator(crossSection);
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(MPS);
  }
  std::map<std::string, std::vector<double> > inits = getParams(MPS);
//auto fs = std::vector<std::function<double(void)> > {};
  std::vector< CorrelatedLL<EventList, pCorrelatedSum&> > csLLs = {};

//auto sigevents_KK = getEvents("signal", pNames, tags[2], dataFile, intFile);
//auto sigMCevents_KK = getEvents("sigMC", pNames, tags[2], dataFile, intFile);
//auto tagevents_KK = getEvents("tag", pNames, tags[2], dataFile, intFile);
//auto tagMCevents_KK = getEvents("tagMC", pNames, tags[2], dataFile, intFile);



//auto events_Kppim = getEvents(pNames, tags[1], dataFile, intFile);
//auto events_Kspipi = getEvents(pNames, tags[2], dataFile, intFile);
//auto types_KK = makeEventTypes(pNames, tags[0]);
//auto types_Kppim = makeEventTypes(pNames, tags[1]);
//auto types_Kspipi = makeEventTypes(pNames, tags[2]);
//auto LL_KK =  tagLL(events_KK, MPS);
//auto LL_Kppim =  tagLL(events_Kppim, MPS);
//auto LL_Kspipi =  tagLL(events_Kspipi, MPS);

//std::stringstream KK_log;
//KK_log<<tags[2]<<"_fit.log";
//std::stringstream KK_fit;
//KK_fit<<tags[2]<<"_plots.root";
//auto KK_logName = KK_log.str();
////auto cs_KK = makeCs(events_KK, MPS);
//auto cs_KK = pCorrelatedSum(sigevents_KK.eventType(), tagevents_KK.eventType(), MPS);
//cs_KK.setEvents(sigevents_KK, tagevents_KK);
//cs_KK.setMC(sigMCevents_KK, tagMCevents_KK);
//cs_KK.prepare();
////auto LL_KK2 = make_likelihood( events_KK["signal"], events_KK["tag"], false, cs_KK);
//auto LL_KK2 = make_likelihood( sigevents_KK, tagevents_KK, false, cs_KK);
//auto mini_KK = Minimiser(LL_KK2, &MPS);
//mini_KK.prepare();
//mini_KK.doFit();
//FitResult * fr_KK = new FitResult(mini_KK);
//FitResult * fr_KK = Fit(LL_KK2, cs_KK, sigevents_KK, tagevents_KK, sigMCevents_KK, tagMCevents_KK, MPS, KK_logName);
//csLLs.push_back(LL_KK2);
INFO("Doing loop of Fits");
for (int i=0; i < tags.size(); i++){
    std::stringstream tag_log;
    auto tagName = split(tags[i],' ')[0];
    tag_log<<tagName<<"_fit.log";
    std::stringstream tag_fit;
    tag_fit<<tagName<<"_plots.root";
    auto tag_plotName = tag_fit.str();
    auto tag_logName = tag_log.str();

    auto sigevents_tag = getEvents("signal", pNames, tags[i], dataFile, intFile);
    auto sigMCevents_tag = getEvents("sigMC", pNames, tags[i], dataFile, intFile);
    auto tagevents_tag = getEvents("tag", pNames, tags[i], dataFile, intFile);
    auto tagMCevents_tag = getEvents("tagMC", pNames, tags[i], dataFile, intFile);

     
    auto cs_tag = pCorrelatedSum(sigevents_tag.eventType(), tagevents_tag.eventType(), MPS);
    cs_tag.setEvents(sigevents_tag, tagevents_tag);
    cs_tag.setMC(sigMCevents_tag, tagMCevents_tag);
    cs_tag.prepare();
    //auto LL_tag2 = make_likelihood( events_tag["signal"], events_tag["tag"], false, cs_tag);
    auto LL_tag2 = make_likelihood( sigevents_tag, tagevents_tag, false, cs_tag);
    auto mini_tag = Minimiser(LL_tag2, &MPS);
    mini_tag.prepare();
    //INFO("Fitting "<<i<<" out of "<<tags.size()  );
    mini_tag.doFit();
  FitResult * fr = new FitResult(mini_tag); 
 int minEvents=15;
auto binning = BinDT( sigevents_tag, MinEvents( minEvents ), Dim( sigevents_tag.eventType().dof() ) );
//void   Chi2Estimator::doChi2( const EventList& sigevents_tag, const EventList& sigMCevents_tag,
//    const std::function<double( const Event& )>& fcn )
//{
  std::vector<Moment> data( binning.size() );
  std::vector<Moment> mc( binning.size() );

  INFO( "Splitting: " << sigevents_tag.size() << " data " << sigMCevents_tag.size() << " amongst " << binning.size()
      << " bins" );
  unsigned int j           = 0;
  double total_data_weight = 0;
  double total_int_weight  = 0;
  for ( auto& d : sigevents_tag ) {
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " data events" );
    double w = d.weight();
    data[binning.getBinNumber( d )].add( d.weight() );
    total_data_weight += w;
    j++;
  }
  j = 0;
  for ( int i=0; i < sigMCevents_tag.size(); i++ ) {
    auto evt1 = sigMCevents_tag[i];
    auto evt2 = tagMCevents_tag[i];
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " sim. events" );
    double w = cs_tag.prob( evt1, evt2 ) * (evt1.weight() / evt1.genPdf()) * (evt2.weight() / evt2.genPdf());
    mc[binning.getBinNumber( evt1 )].add( w );
    total_int_weight += w;
    j++;
  }
  double chi2 = 0;

  for ( unsigned int i = 0; i < binning.size(); ++i ) {
    mc[i].rescale( total_data_weight / total_int_weight );
    double delta = data[i].val() - mc[i].val();
    double tChi2 = delta * delta / ( data[i].val() + mc[i].var() );
    chi2 += tChi2;
  }

  auto Bins = binning.size();
  auto fitFractionss = cs_tag.fitFractions( fr->getErrorPropagator() ); 
  std::vector<FitFraction> ffs;
  for (auto fitFractions : fitFractionss){
    for (auto fitFraction : fitFractions){
      ffs.push_back(fitFraction);
    }
  }
  fr->addFractions( ffs );
  //chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2, Bins );
  fr->print();
  
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(tag_logName, pulls);


   if (doProjections){
    INFO( "norm[1] = " << cs_tag.norm() );
    TFile* f = TFile::Open(tag_plotName.c_str(),"RECREATE");

    auto projections = sigevents_tag.eventType().defaultProjections(nBins);
    for( auto& projection : projections ){
      auto data_plot = projection(sigevents_tag);
      auto hist = projection.plot();
      for(unsigned i = 0 ; i != sigMCevents_tag.size(); ++i)
      {
        hist->Fill( projection( sigMCevents_tag[i] ), cs_tag.prob( sigMCevents_tag[i], tagMCevents_tag[i] ) );
      }
      hist->Scale( data_plot->Integral() / hist->Integral() );
      hist->SetName( (std::string("MC_")+hist->GetName()).c_str() );
      hist->Write();
      data_plot->Write();
      auto pull = (TH1D*) hist->Clone();
      pull->Add(data_plot, -1);
      pull->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
      pull->Write();
    }
    auto p2 = sigMCevents_tag.eventType().defaultProjections(nBins);
    for( unsigned i = 0 ; i != p2.size() -1; ++i )
    {
      for( unsigned j=i+1; j < p2.size(); ++j )
      {
        auto dalitz = Projection2D( p2[i], p2[j] );
        auto hdalitz = dalitz.plot();
        auto data_plot = sigevents_tag.makeProjection(dalitz);
        for( unsigned event = 0 ; event != sigMCevents_tag.size(); ++event )
        {
          auto pos = dalitz(sigMCevents_tag[event]);
          hdalitz->Fill( pos.first, pos.second, cs_tag.prob( sigMCevents_tag[event], tagMCevents_tag[event] ) );
        }
        hdalitz->Scale( data_plot->Integral() / hdalitz->Integral() );
        hdalitz->SetName( ( std::string("MC_") + hdalitz->GetName() ).c_str() );
        hdalitz->Write();

        data_plot->Write(); 
        auto pull2D = (TH2D*) hdalitz->Clone();
        pull2D->Add(data_plot, -1);

        pull2D->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
        pull2D->Write();
      }
    }

    f->Close();
    
    }

   fr->writeToFile(tag_logName);
    
    //FitResult * fr_tag = Fit(LL_tag2, cs_tag, sigevents_tag, tagevents_tag, sigMCevents_tag, tagMCevents_tag, MPS, tag_logName, inits);
    csLLs.push_back(LL_tag2);

}



if (doCombFit){
//doPlots(KK_fit.str(), cs_KK, events_KK, 100);

//auto LLs = {LL_KK2, LL_Kppim2};
//INFO("Trying to do combined fit");
// auto combLL = SumLL<CorrelatedLL<EventList, pCorrelatedSum&> >({csLLs[0], csLLs[1]});
// auto combMini = Minimiser(combLL, &MPS);
// INFO("Minimising now");
// combMini.doFit();
//FitResult * fr_KK = new FitResult(mini_KK);
  std::vector<CorrelatedLL<EventList, pCorrelatedSum&> > LLs;
//FitResult * fr_comb = new FitResult(combMini);
  if (tags.size()>0){
    INFO("Making events 0");
    auto sigevents_0 = getEvents("signal", pNames, tags[0], dataFile, intFile);
    auto sigMCevents_0 = getEvents("sigMC", pNames, tags[0], dataFile, intFile);
    auto tagevents_0 = getEvents("tag", pNames, tags[0], dataFile, intFile);
    auto tagMCevents_0 = getEvents("tagMC", pNames, tags[0], dataFile, intFile);
    INFO("Making pCorrelatedSum 0");
    auto cs_0 = pCorrelatedSum(sigevents_0.eventType(), tagevents_0.eventType(), MPS);
    cs_0.setEvents(sigevents_0, tagevents_0);
    cs_0.setMC(sigMCevents_0, tagMCevents_0);
    cs_0.prepare();
    INFO("Making LL 0");
    //auto LL_tag2 = make_likelihood( events_tag["signal"], events_tag["tag"], false, cs_tag);
    auto LL_0 = make_likelihood( sigevents_0, tagevents_0, false, cs_0);
    auto mini_0 = Minimiser(LL_0, &MPS);
    LLs = {LL_0};
    

    if (tags.size() > 1){
      INFO("Making events 1");
      auto sigevents_1 = getEvents("signal", pNames, tags[1], dataFile, intFile);
      auto sigMCevents_1 = getEvents("sigMC", pNames, tags[1], dataFile, intFile);
      auto tagevents_1 = getEvents("tag", pNames, tags[1], dataFile, intFile);
      auto tagMCevents_1 = getEvents("tagMC", pNames, tags[1], dataFile, intFile);
      INFO("Making pCorrelatedSum 1");
      auto cs_1 = pCorrelatedSum(sigevents_1.eventType(), tagevents_1.eventType(), MPS);
      cs_1.setEvents(sigevents_1, tagevents_1);
      cs_1.setMC(sigMCevents_1, tagMCevents_1);
      cs_1.prepare();
      INFO("Making LL 1");
      auto LL_1 = make_likelihood( sigevents_1, tagevents_1, false, cs_1);
      auto mini_1 = Minimiser(LL_1, &MPS);
      LLs = {LL_0, LL_1};
     
      if (tags.size() > 2){
        INFO("Making events 2");
        auto sigevents_2 = getEvents("signal", pNames, tags[2], dataFile, intFile);
        auto sigMCevents_2 = getEvents("sigMC", pNames, tags[2], dataFile, intFile);
        auto tagevents_2 = getEvents("tag", pNames, tags[2], dataFile, intFile);
        auto tagMCevents_2 = getEvents("tagMC", pNames, tags[2], dataFile, intFile);
        INFO("Making pCorrelatedSum 2");
        auto cs_2 = pCorrelatedSum(sigevents_2.eventType(), tagevents_2.eventType(), MPS);
        cs_2.setEvents(sigevents_2, tagevents_2);
        cs_2.setMC(sigMCevents_2, tagMCevents_2);
        cs_2.prepare(); 
        INFO("Making LL 2");
        auto LL_2 = make_likelihood( sigevents_2, tagevents_2, false, cs_2);
        auto mini_2 = Minimiser(LL_2, &MPS);
        LLs = {LL_0, LL_1, LL_2};
        }
      }
    }
  else{
        INFO("No tags - returning 0");
        return 0;
      }






  INFO("Making combLL"); 
    auto combLL = SumLL<CorrelatedLL<EventList, pCorrelatedSum&> >(LLs);
    auto combMini = Minimiser(combLL, &MPS);
    INFO("Minimising now");
    combMini.doFit();

  FitResult * fr = new FitResult(combMini); 
  fr->print();
   fr->writeToFile(logFile);
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(logFile, pulls);




   fr->writeToFile(logFile);
 }

/*
 for( auto& tag : tags )
 {
    EventType signalType( pNames );
    auto tokens       = split(tag, ' ');
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType(); 
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");

    EventList sigEvents;
    EventList tagEvents;
    EventList sigMCEvents;
    EventList tagMCEvents;

    if (QcGen2){
        INFO("Generated Events used QcGen2");
    std::stringstream signame;
    signame<<":Signal_";
    signame<<tokens[0];
    std::stringstream tagname;
    tagname<<":Tag_";
    tagname<<tokens[0];
      sigEvents = EventList(dataFile + signame.str() , signalType);
      tagEvents = EventList(dataFile + tagname.str() , tagType);
      sigMCEvents = EventList(intFile + signame.str() , signalType);
      tagMCEvents = EventList(intFile + tagname.str() , tagType);

    }
    else{

        INFO("Generated Events used QcGenerator");
    sigEvents = EventList(dataFile +":"+ tokens[0], signalType, Branches(sigBranches));
    sigMCEvents = EventList(intFile +":"+ tokens[0], signalType, Branches(sigBranches));
    tagEvents = EventList(dataFile +":"+ tokens[0], tagType, Branches(tagBranches));
    tagMCEvents = EventList(intFile +":"+ tokens[0], tagType, Branches(tagBranches));
    }
    for( auto& event : sigEvents ) event.setGenPdf(1) ;
    for( auto& event : tagEvents ) event.setGenPdf(1) ;
    
  //   CorrelatedSum cs(signalType, tagType, MPS);
  //  cs.setEvents(sigEvents, tagEvents);
  //  cs.setMC(sigMCEvents, tagMCEvents);
  //   cs.prepare();  
  //  auto csLL = make_likelihood(sigEvents, tagEvents, cs);
  //  csLL.setEvents(sigEvents, tagEvents);
    

    //INFO( "norm[0] = " << cs.norm() );
   
  //  cs.debugNorm();
// Minimiser mini( csLL, &MPS );
   // for (int i=0;  i<nFits; i++){
   
    
    mini.GradientTest();
   // for (int i=0 ; i < nFits; i++){
    //   mini.setTolerance(std::pow(10, 2 - i));
//   if (doFit) mini.doFit();
   // }
   // 

   
    
    pCorrelatedSum cs(signalType, tagType, MPS);
    cs.setEvents(sigEvents, tagEvents);
    cs.setMC(sigMCEvents, tagMCEvents);
     cs.prepare();  
//    auto csLL = make_likelihood(sigEvents, tagEvents, false, cs);
    CorrelatedLL<EventList, pCorrelatedSum&> csLL =  make_likelihood(sigEvents, tagEvents, false, cs);
    csLL.setEvents(sigEvents, tagEvents);
    INFO( "norm[0] = " << cs.norm() );
    cs.debugNorm();

    CorrelatedLL<EventList, pCorrelatedSum&> * csLL1 = new CorrelatedLL<EventList, pCorrelatedSum&>(false, cs);
   // csLLs.push_back(csLL);
    //csLLs.push_back(csLL);
    //

  std::vector< CorrelatedLL<EventList, pCorrelatedSum&> > LLs = {csLL, csLL};

 //   auto pCsLL=*csLL;
 Minimiser mini( csLL, &MPS );

auto sumLLs= SumLL< CorrelatedLL<EventList, pCorrelatedSum&> >(LLs);

 auto combMini = Minimiser(sumLLs, &MPS);
 combMini.prepare();
// combMini.doFit();

 INFO("type of CS = "<<typeid(cs).name());
 INFO("type of CSLL = "<<typeid(csLLs).name());
 INFO("CLL = "<<csLL.getVal());
   // for (int i=0;  i<nFits; i++){
//fs.push_back(mini.fitFunction());
   
    
    mini.GradientTest();
   // for (int i=0 ; i < nFits; i++){
    //   mini.setTolerance(std::pow(10, 2 - i));
   //if (doFit) mini.doFit();
   // }
   // 
      
     

//    }
    FitResult * fr = new FitResult(mini);
 auto dataEvents = sigEvents;
 auto mcEvents = sigMCEvents;
 int minEvents=15;
auto binning = BinDT( dataEvents, MinEvents( minEvents ), Dim( dataEvents.eventType().dof() ) );
//void   Chi2Estimator::doChi2( const EventList& dataEvents, const EventList& mcEvents,
//    const std::function<double( const Event& )>& fcn )
//{
  std::vector<Moment> data( binning.size() );
  std::vector<Moment> mc( binning.size() );

  INFO( "Splitting: " << dataEvents.size() << " data " << mcEvents.size() << " amongst " << binning.size()
      << " bins" );
  unsigned int j           = 0;
  double total_data_weight = 0;
  double total_int_weight  = 0;
  for ( auto& d : dataEvents ) {
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " data events" );
    double w = d.weight();
    data[binning.getBinNumber( d )].add( d.weight() );
    total_data_weight += w;
    j++;
  }
  j = 0;
  for ( int i=0; i < mcEvents.size(); i++ ) {
    auto evt1 = mcEvents[i];
    auto evt2 = tagMCEvents[i];
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " sim. events" );
    double w = cs.prob( evt1, evt2 ) * (evt1.weight() / evt1.genPdf()) * (evt2.weight() / evt2.genPdf());
    mc[binning.getBinNumber( evt1 )].add( w );
    total_int_weight += w;
    j++;
  }
  double chi2 = 0;

  for ( unsigned int i = 0; i < binning.size(); ++i ) {
    mc[i].rescale( total_data_weight / total_int_weight );
    double delta = data[i].val() - mc[i].val();
    double tChi2 = delta * delta / ( data[i].val() + mc[i].var() );
    chi2 += tChi2;
  }
  
  auto Bins = binning.size();
  auto fitFractionss = cs.fitFractions( fr->getErrorPropagator() ); 
  std::vector<FitFraction> ffs;
  for (auto fitFractions : fitFractionss){
    for (auto fitFraction : fitFractions){
      ffs.push_back(fitFraction);
    }
  }
  fr->addFractions( ffs );
  //chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2, Bins );
  fr->print();
   fr->writeToFile(logFile);
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(logFile, pulls);


    if (doProjections){
    INFO( "norm[1] = " << cs.norm() );
    TFile* f = TFile::Open(plotFile.c_str(),"RECREATE");

    auto projections = signalType.defaultProjections(nBins);
    for( auto& projection : projections ){
      auto data_plot = projection(sigEvents);
      auto hist = projection.plot();
      for(unsigned i = 0 ; i != sigMCEvents.size(); ++i)
      {
        hist->Fill( projection( sigMCEvents[i] ), cs.prob( sigMCEvents[i], tagMCEvents[i] ) );
      }
      hist->Scale( data_plot->Integral() / hist->Integral() );
      hist->SetName( (std::string("MC_")+hist->GetName()).c_str() );
      hist->Write();
      data_plot->Write();
      auto pull = (TH1D*) hist->Clone();
      pull->Add(data_plot, -1);
      pull->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
      pull->Write();
    }
    auto p2 = signalType.defaultProjections(nBins);
    for( unsigned i = 0 ; i != p2.size() -1; ++i )
    {
      for( unsigned j=i+1; j < p2.size(); ++j )
      {
        auto dalitz = Projection2D( p2[i], p2[j] );
        auto hdalitz = dalitz.plot();
        auto data_plot = sigEvents.makeProjection(dalitz);
        for( unsigned event = 0 ; event != sigMCEvents.size(); ++event )
        {
          auto pos = dalitz(sigMCEvents[event]);
          hdalitz->Fill( pos.first, pos.second, cs.prob( sigMCEvents[event], tagMCEvents[event] ) );
        }
        hdalitz->Scale( data_plot->Integral() / hdalitz->Integral() );
        hdalitz->SetName( ( std::string("MC_") + hdalitz->GetName() ).c_str() );
        hdalitz->Write();

        data_plot->Write(); 
        auto pull2D = (TH2D*) hdalitz->Clone();
        pull2D->Add(data_plot, -1);

        pull2D->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
        pull2D->Write();
      }
    }

    f->Close();
    
    }
    
  }
 //auto combLL = SumLL<CorrelatedLL<EventList, pCorrelatedSum&> >(csLLs);
 
 //auto combLL = SumLL<std::function<double(void)> >(fs);
// auto combMini = Minimiser(combLL, &MPS);
 //combMini.prepare();
//auto ll0 = combLL.getVal();
//INFO("LL = "<<ll0);
 //combMini.doFit();
 //INFO("Number of LL = "<<csLLs.size());
// double x = csLLs[0]->getVal();
// INFO("x = "<<x);
  //for (auto & csll : csLLs){
//      auto l = csll.getVal();
//      INFO("l = "<<l);
  //}
  */
  return 0;
}


std::vector<std::string> makeBranches(EventType Type, std::string prefix){
  auto n = Type.finalStates().size();
  std::vector<std::string> branches;
  std::vector<std::string> varNames = {"PX", "PY", "PZ", "E"};
  for (long unsigned int i=0; i<n; i++){
    auto part = replaceAll(Type.finalStates()[i], "+", "p");
    part = replaceAll(part, "-", "m");
    for (auto varName : varNames){
      std::ostringstream stringStream;
      stringStream<<prefix<<part<<"_"<<varName;
      branches.push_back(stringStream.str());
    }
  }
  return branches;
}

template <typename PDF>
FitResult* Fit( PDF&& ll, pCorrelatedSum& cs, EventList data1, EventList data2, EventList mc1, EventList mc2, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  MPS.loadFromStream();

  ll.setEvents( data1, data2 );
 
//  pdf.setMC(mc1, mc2);
  /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
  Minimiser mini( ll, &MPS );
   
  mini.prepare();
  mini.GradientTest();
 INFO("Fitting now"); 
  mini.doFit();
  INFO("Making Fit Result");
  FitResult* fr = new FitResult(mini);

//  std::map<std::string, std::vector<double> > inits = getParams(MPS);
  // Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
  //   The structure assumed the PDF is some SumPDF<T1,T2,...>. 
 

 auto dataEvents = data1;
 auto mcEvents = mc1;
 auto tagMCEvents = mc2;


 int minEvents=15;
auto binning = BinDT( dataEvents, MinEvents( minEvents ), Dim( dataEvents.eventType().dof() ) );
//void   Chi2Estimator::doChi2( const EventList& dataEvents, const EventList& mcEvents,
//    const std::function<double( const Event& )>& fcn )
//{
  std::vector<Moment> data( binning.size() );
  std::vector<Moment> mc( binning.size() );

  INFO( "Splitting: " << dataEvents.size() << " data " << mcEvents.size() << " amongst " << binning.size()
      << " bins" );
  unsigned int j           = 0;
  double total_data_weight = 0;
  double total_int_weight  = 0;
  for ( auto& d : dataEvents ) {
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " data events" );
    double w = d.weight();
    data[binning.getBinNumber( d )].add( d.weight() );
    total_data_weight += w;
    j++;
  }
  j = 0;
  for ( int i=0; i < mcEvents.size(); i++ ) {
    auto evt1 = mcEvents[i];
    auto evt2 = tagMCEvents[i];
    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " sim. events" );
    double w = cs.prob( evt1, evt2 ) * (evt1.weight() / evt1.genPdf()) * (evt2.weight() / evt2.genPdf());
    mc[binning.getBinNumber( evt1 )].add( w );
    total_int_weight += w;
    j++;
  }
  double chi2 = 0;

  for ( unsigned int i = 0; i < binning.size(); ++i ) {
    mc[i].rescale( total_data_weight / total_int_weight );
    double delta = data[i].val() - mc[i].val();
    double tChi2 = delta * delta / ( data[i].val() + mc[i].var() );
    chi2 += tChi2;
  }

  auto Bins = binning.size();
  auto fitFractionss = cs.fitFractions( fr->getErrorPropagator() ); 
  std::vector<FitFraction> ffs;
  for (auto fitFractions : fitFractionss){
    for (auto fitFraction : fitFractions){
      ffs.push_back(fitFraction);
    }
  }
  fr->addFractions( ffs );
  //chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2, Bins );
  fr->print();
   fr->writeToFile(logFile);
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(logFile, pulls);



/* 
  for_each(pdf.m_amps, [&]( auto& f ){
    std::function<double(const Event&, const Event&)> FCN_sig = [&](const Event& event1, const Event& event2){ return f.prob_unnormalised(event1, event2) ; };
    auto mc_plot3 = mc1.makeDefaultProjections(WeightFunction(f), Prefix("Model_cat"+std::to_string(counter)));
    for( auto& plot : mc_plot3 )
    {
      plot->Scale( ( data1.integral()) / plot->Integral() );
      plot->Write();
    }
    counter++;
  } );
  */
  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
 /* 
  Chi2Estimator chi2( data1, mc1, pdf, 15 );
  chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2.chi2(), chi2.nBins() );
  */
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );

  fr->print();
  return fr;
}

void add_CP_conjugate( MinuitParameterSet& mps )
{
  std::vector<MinuitParameter*> tmp;
  for( auto& param : mps ){
    const std::string name = param->name();
    size_t pos=0;
    std::string new_name = name; 
    int sgn=1;
    Flag flag = param->flag();
    //Flag flag = Flag::Fix;
    if( name.find("::") != std::string::npos ){
      pos = name.find("::");
      auto props = AmpGen::ParticlePropertiesList::get( name.substr(0,pos), true );
      if( props != 0 ) new_name = props->anti().name() + name.substr(pos); 
    }
    else { 
      auto tokens=split(name,'_');
      std::string reOrIm = *tokens.rbegin();
      std::string name   = tokens[0];
      if ( reOrIm == "Re" || reOrIm == "Im" ){
        auto p = Particle( name ).conj();
        sgn = reOrIm == "Re" ? p.quasiCP() : 1; 
        new_name = p.uniqueString() +"_"+reOrIm;
      }
      else if( tokens.size() == 2 ) {
        auto props = AmpGen::ParticlePropertiesList::get( name );
        if( props != 0  ) new_name = props->anti().name() + "_" + tokens[1]; 
      }
    }
    if( mps.find( new_name ) == nullptr ){
      tmp.push_back( new MinuitParameter(new_name, flag, sgn * param->mean(), param->err(), 0, 0));
    }
  }
  for( auto& p : tmp ) mps.add( p );
}

std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps){
  std::map<std::string, std::vector<double> > out;
  for (auto& param : mps){
    std::vector<double> vect = {param->mean(), param->err()};
    out[param->name()] = vect;
  }
  return out;
}


MinuitParameterSet * copyMPS(MinuitParameterSet& mps){
  MinuitParameterSet * newMPS = new MinuitParameterSet; 
  for (auto& param : mps){
    newMPS->add(param);
  }
  return newMPS;
}

std::map<std::string, double > getPulls(std::map<std::string, std::vector<double> > fits, std::map<std::string, std::vector<double> > inits)  {
  std::string output = "";
 std::map<std::string, double> out;
 for(map<std::string, std::vector<double> >::iterator init = inits.begin(); init != inits.end(); ++init) {
  for(map<std::string, std::vector<double> >::iterator fit = fits.begin(); fit != fits.end(); ++fit) {
    std::string initName = init->first;
    std::vector<double> initParams = init->second;
    std::string fitName = fit->first;
    std::vector<double> fitParams = fit->second;
    if (initName == fitName){
        double pull = fitParams[0] - initParams[0];
        if (fitParams[1] != 0){
          pull /= pow(pow(fitParams[1], 2) + pow(initParams[1],2),0.5);
        }
        out[fitName] = pull;
    }
  }
 }
  return out;

}

void writePulls(std::string fileName, std::map<std::string, double> pulls){
  std::ofstream outfile;
  outfile.open(fileName, std::ios_base::app);
      for(map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
        outfile<<"Pull "<<it->first<<" "<<it->second<<"\n";
       }

}



template <typename pdftype>
void doPlots(std::string plotFile, pdftype&& cs,  EventList sigEvents, EventList tagEvents, EventList sigMCEvents, EventList tagMCEvents , int nBins)
{
    

    auto signalType = sigEvents.eventType();
    auto tagType = tagEvents.eventType();
    INFO( "norm[1] = " << cs.norm() );
    TFile* f = TFile::Open(plotFile.c_str(),"RECREATE");

    auto projections = signalType.defaultProjections(nBins);
    for( auto& projection : projections ){
      auto data_plot = projection(sigEvents);
      auto hist = projection.plot();
      for(unsigned i = 0 ; i != sigMCEvents.size(); ++i)
      {
        hist->Fill( projection( sigMCEvents[i] ), cs.prob( sigMCEvents[i], tagMCEvents[i] ) );
      }
      hist->Scale( data_plot->Integral() / hist->Integral() );
      hist->SetName( (std::string("MC_")+hist->GetName()).c_str() );
      hist->Write();
      data_plot->Write();
      auto pull = (TH1D*) hist->Clone();
      pull->Add(data_plot, -1);
      pull->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
      pull->Write();
    }
    auto p2 = signalType.defaultProjections(nBins);
    for( unsigned i = 0 ; i != p2.size() -1; ++i )
    {
      for( unsigned j=i+1; j < p2.size(); ++j )
      {
        auto dalitz = Projection2D( p2[i], p2[j] );
        auto hdalitz = dalitz.plot();
        auto data_plot = sigEvents.makeProjection(dalitz);
        for( unsigned event = 0 ; event != sigMCEvents.size(); ++event )
        {
          auto pos = dalitz(sigMCEvents[event]);
          hdalitz->Fill( pos.first, pos.second, cs.prob( sigMCEvents[event], tagMCEvents[event] ) );
        }
        hdalitz->Scale( data_plot->Integral() / hdalitz->Integral() );
        hdalitz->SetName( ( std::string("MC_") + hdalitz->GetName() ).c_str() );
        hdalitz->Write();

        data_plot->Write(); 
        auto pull2D = (TH2D*) hdalitz->Clone();
        pull2D->Add(data_plot, -1);

        pull2D->SetName( (std::string("Pull_") + data_plot->GetName()).c_str() );
        pull2D->Write();
      }
    }

    f->Close();
}


///CorrelatedLL
CorrelatedLL<EventList, pCorrelatedSum&> tagLL(std::map<std::string, EventList> events, MinuitParameterSet& MPS){

INFO("Test function");

    EventList sigEvents = events["signal"];
    EventList tagEvents = events["tag"];
    EventList sigMCEvents = events["sigMC"];
    EventList tagMCEvents = events["tagMC"];

    auto signalType = sigEvents.eventType();
    auto tagType = tagEvents.eventType();


    for( auto& event : sigEvents ) event.setGenPdf(1) ;
    for( auto& event : tagEvents ) event.setGenPdf(1) ;

   
    
    pCorrelatedSum cs(signalType, tagType, MPS);
    cs.setEvents(sigEvents, tagEvents);
    cs.setMC(sigMCEvents, tagMCEvents);
    cs.prepare();  

    CorrelatedLL<EventList, pCorrelatedSum&> csLL =  make_likelihood(sigEvents, tagEvents, false, cs);
    csLL.setEvents(sigEvents, tagEvents);
    INFO( "norm[0] = " << cs.norm() );
    cs.debugNorm();

    return csLL;
}


EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile, std::string intFile){
   
    /*
    EventType signalType( sigName );
    INFO("Getting Tag name split "<<tagName);
    auto tokens       = split(tagName, ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto tagParticle  = Particle(tokens[1], {}, false);
    INFO("Build EventType");
    EventType    tagType = tagParticle.eventType(); 
    */

    INFO("Getting Tag name split "<<tagName);
    auto tokens       = split(tagName, ' ');
    INFO("Building Tag particle"<<tokens[1]);
    auto types = makeEventTypes(sigName, tagName);
    INFO("Build EventType");
    auto signalType = types["signal"];
    auto tagType = types["tag"];
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");

    EventList sigEvents;
    EventList tagEvents;
    EventList sigMCEvents;
    EventList tagMCEvents;
    EventList null;


    INFO("Generated Events used QcGen2");
    std::stringstream signame;
    signame<<":Signal_";
    signame<<tokens[0];
    std::stringstream tagname;
    tagname<<":Tag_";
    tagname<<tokens[0];
    sigEvents = EventList(dataFile + signame.str() , signalType);
    tagEvents = EventList(dataFile + tagname.str() , tagType);
    sigMCEvents = EventList(intFile + signame.str() , signalType);
    tagMCEvents = EventList(intFile + tagname.str() , tagType);
    
    if (type=="signal"){
     return sigEvents;
    }
    else if (type=="tag"){ 
        return tagEvents;
    }
    else if (type=="sigMC"){
        return sigMCEvents;
    }
    else if (type=="tagMC") {
        return tagMCEvents;
    }
    else {
        return null;
    }
    


}

std::map<std::string, EventType> makeEventTypes(std::vector<std::string> sigName, std::string tagName){
  EventType signalType( sigName );
  INFO("Getting Tag name split "<<tagName);
  auto tokens       = split(tagName, ' ');
  INFO("Building Tag particle"<<tokens[1]);
  auto tagParticle  = Particle(tokens[1], {}, false);
  INFO("Build EventType");
  EventType    tagType = tagParticle.eventType(); 
  auto sigBranches = makeBranches(signalType, "");
  auto tagBranches = makeBranches(tagType, "Tag_");

  std::map<std::string, EventType> types = {};

  types["signal"] = signalType;
  types["tag"] = tagType;
  return types;

}

pCorrelatedSum makeCs(std::map<std::string, EventList> events,  MinuitParameterSet& MPS){
  

    auto sigEvents = events["signal"];
    auto sigMCEvents = events["sigMC"];
    auto tagEvents = events["tag"];
    auto tagMCEvents = events["tagMC"];
    auto signalType = sigEvents.eventType();
    auto tagType = tagEvents.eventType();
    pCorrelatedSum cs(signalType, tagType, MPS);
    cs.setEvents(sigEvents, tagEvents);
    cs.setMC(sigMCEvents, tagMCEvents);
    cs.prepare();  
    return cs;
   
}
