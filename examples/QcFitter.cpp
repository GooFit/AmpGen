#include "AmpGen/Psi3770.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/corrEventList.h"
//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

std::vector<std::string> makeBranches(EventType Type, std::string prefix);
template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data1, EventList& data2, EventList& mc1, EventList& mc2, MinuitParameterSet& MPS );

void add_CP_conjugate( MinuitParameterSet& mps );
template <typename pdftype>
void doPlots(pdftype&& pdf, EventList& data, EventList& mc, std::string plotFile);
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
    int nFits = NamedParameter<int>("nFits", 1, "number of repeats of mini.doFits() for debug purposes!");


  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
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
  auto yc = DTYieldCalculator(crossSection);
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(MPS);
  }
 for( auto& tag : tags )
 {
    EventType signalType( pNames );
    auto tokens       = split(tag, ' ');
    auto tagParticle  = Particle(tokens[1], {}, false);
    EventType    tagType = tagParticle.eventType(); 
    auto sigBranches = makeBranches(signalType, "");
    auto tagBranches = makeBranches(tagType, "Tag_");
    EventList sigEvents(dataFile +":"+ tokens[0], signalType, Branches(sigBranches));
    EventList sigMCEvents(intFile +":"+ tokens[0], signalType, Branches(sigBranches));
    EventList tagEvents(dataFile +":"+ tokens[0], signalType, Branches(tagBranches));
    EventList tagMCEvents(intFile +":"+ tokens[0], signalType, Branches(tagBranches));
    for( auto& event : sigEvents ) event.setGenPdf(1) ;
    for( auto& event : tagEvents ) event.setGenPdf(1) ;
    
    CorrelatedSum cs(signalType, tagType, MPS);
    cs.setEvents(sigEvents, tagEvents);
    cs.setMC(sigMCEvents, tagMCEvents);
    auto csLL = make_likelihood(sigEvents, tagEvents, cs);
    csLL.setEvents(sigEvents, tagEvents);
    
    cs.prepare(); 
    INFO( "norm[0] = " << cs.norm() );

    Minimiser mini( csLL, &MPS );
    
    for (int i=0 ; i < nFits; i++) mini.doFit();
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

  //chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2, Bins );
  fr->print();
   fr->writeToFile(logFile);



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
      }
    }

    f->Close();
  }
  
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
FitResult* doFit( PDF&& ll, EventList& data1, EventList& mc1, EventList& data2, EventList& mc2, MinuitParameterSet& MPS )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();

  ll.setEvents( data1, data2 );
 
//  pdf.setMC(mc1, mc2);
  /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
  Minimiser mini( ll, &MPS );
   
  mini.prepare();
  mini.GradientTest();
  
  mini.doFit();
  FitResult* fr = new FitResult(mini);

  // Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
  //   The structure assumed the PDF is some SumPDF<T1,T2,...>. 
 

  corrEventList cEL_data(data1, data2);
  corrEventList cEL_mc(mc1, mc2);
  unsigned int counter= 1;
  for_each(ll.pdfs(), [&](auto& pdf){
   INFO("Printing all the pdfs in my likelihood");
    auto pfx1 = Prefix("Model_cat"+std::to_string(counter));
    //
    //For our QcFitter we can't just ask for projections like in SignalOnlyFitter, 
    //we need to have two sets of projections for each event.
    //
    ////std::vector<std::vector<TH1D*> > hists_data =  cEL_data.makeDefaultProjections();
    std::vector<std::vector<TH1D*> > hists_mc =  cEL_mc.makeDefaultProjections();
    //auto mc1_plot3 = mc1.makeDefaultProjections(WeightFunction(pdf), pfx1);
    counter++;
    
  });

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
      tmp.push_back( new MinuitParameter(new_name, Flag::Free, sgn * param->mean(), param->err(), 0, 0));
    }
  }
  for( auto& p : tmp ) mps.add( p );
}

template <typename pdftype>
void doPlots(pdftype&& pdf, EventList& data, EventList& mc, std::string plotFile){
  TFile * output = TFile::Open(plotFile.c_str(), "RECREATE");
  output->cd();

   /*
  for (auto pdf : ll.pdfs()){
    auto pfx = Prefix("Model_cat"+std::to_string(counter));
    auto mc_plot3 = mc.makeDefaultProjections(WeightFunction(pdf), pfx);
    for (auto& plot : mc_plot3){
      plot->Scale( (data.integral())/plot->Integral() );
      plot->Write();
    }
    counter++;

  }
  */
  output->Close();

}
