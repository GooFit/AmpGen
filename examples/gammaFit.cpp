#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/CombCorrLL.h"
#include "AmpGen/CombLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>

//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

std::vector<std::string> makeBranches(EventType Type, std::string prefix);
//template <typename PDF>
//FitResult* Fit( PDF&& pdf, pCorrelatedSum& cs, std::map<std::string, EventList> events, MinuitParameterSet& MPS , std::string logFile);
template <typename PDF>
FitResult* Fit( PDF&& ll, pCorrelatedSum& cs, EventList data1, EventList data2, EventList mc1, EventList mc2, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits, int maxAttempts , bool repeatFits);

template <typename PDF>
FitResult* Fit( PDF&& ll, pCoherentSum& cs, EventList data1, EventList mc1, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits, int maxAttempts, bool repeatFits );

void FitTag(EventList& sigevents_tag, EventList& tagevents_tag, EventList& sigMCevents_tag, EventList& tagMCevents_tag, MinuitParameterSet* MPS_tag, int maxAttempts, std::map<std::string,  std::vector<double> > inits  , std::string tag_plotName, std::string tag_logName, int nBins, bool doProjections);


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
void doPlots(std::string plotFile, pdftype&& cs,  EventList sigEvents, EventList tagEvents, EventList sigMCEvents, EventList tagMCEvents , int nBins);


template <typename pdftype>
void doPlots(std::string plotFile, pdftype&& cs,  EventList sigEvents,  EventList sigMCEvents, int nBins);

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
  bool doCombFit      = NamedParameter<bool>("doCombFit", true, "Do combined fit");
  bool doTagFit      = NamedParameter<bool>("doTagFit", true, "Do fit for each tag");
  int  maxAttempts = NamedParameter<int>("maxAttempts", 5, "Max attempts to get a valid minimum from Minuit2");
  bool QcGen2 = NamedParameter<bool>("QcGen2", false, "internal boolean - for new QcGenerator");
  bool doFit = NamedParameter<bool>("doFit", true, "Do the fit");
  bool repeatFits = NamedParameter<bool>("repeatFits", false);
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

INFO("Doing loop of Fits");

 std::vector<EventList> SigData;
 std::vector<EventList> TagData;
 std::vector<EventList> SigInt;
 std::vector<EventList> TagInt;
 std::vector<EventType> SigType;
 std::vector<EventType> TagType;
 std::vector<std::string> sfList;

 std::vector<CorrelatedLL<EventList, pCorrelatedSum&> > totalLL;
// SimFit totalLL;


  bool doCorrFit = NamedParameter<bool>("doCorrFit", false);
 if (doCorrFit){
for (int i=0; i < tags.size(); i++){
MinuitParameterSet * MPS_tag = new MinuitParameterSet();
  MPS_tag->loadFromStream();
  if (makeCPConj){
    INFO("Making CP conjugate states");
    add_CP_conjugate(*MPS_tag);
  }
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

    SigData.push_back(sigevents_tag);
    TagData.push_back(tagevents_tag);
    SigInt.push_back(sigMCevents_tag);
    TagInt.push_back(tagMCevents_tag);
    SigType.push_back(sigevents_tag.eventType());
    TagType.push_back(tagevents_tag.eventType());
    sfList.push_back("Psi3770");

      auto cs_tag = pCorrelatedSum(sigevents_tag.eventType(), tagevents_tag.eventType(), *MPS_tag, "Psi3770");
      cs_tag.setEvents(sigevents_tag, tagevents_tag);
      cs_tag.setMC(sigMCevents_tag, tagMCevents_tag);
      cs_tag.prepare();
      //auto LL_tag2 = make_likelihood( events_tag["signal"], events_tag["tag"], false, cs_tag);
      auto LL_tag2 = make_likelihood( sigevents_tag, tagevents_tag,cs_tag);
      totalLL.push_back(LL_tag2);
    if (doTagFit) { 
      FitResult * fr = Fit(LL_tag2, cs_tag, sigevents_tag, tagevents_tag, sigMCevents_tag, tagMCevents_tag, MPS, tag_logName, inits, maxAttempts, repeatFits);
      doPlots(tag_plotName, cs_tag, sigevents_tag, tagevents_tag, sigMCevents_tag, tagMCevents_tag, nBins);
      fr->writeToFile(tag_logName);
    }     
  }



  if (doCombFit){
    CombCorrLL combLL = CombCorrLL(SigData, TagData, SigInt, TagInt, SigType, TagType, MPS, sfList);
    auto combLL2 = SumLL<CorrelatedLL<EventList, pCorrelatedSum&>>(totalLL);
    INFO("CombCorrLL = "<<combLL.getVal());
//    auto commLL2 = SumLL(_LLs);
    INFO("Making Combined Minimiser object");
//    INFO("totalLL = "<<totalLL.getVal());
    Minimiser combMini = Minimiser(combLL, &MPS);
    combMini.gradientTest();
   // combMini.prepare();
    INFO("Minimising now");
    int attempt = 1;
      combMini.doFit(); 
      if (repeatFits){
      if (combMini.status() != 0){
        INFO("Didn't seem to get a minimum (returned "<<combMini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
      while (attempt < maxAttempts && combMini.status() != 0){
        INFO("Didn't seem to get a minimum (returned "<<combMini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
        combMini.doFit();
        attempt++;
      }
      }
      }

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
 }

  auto Btags           = NamedParameter<std::string>("BTagTypes" , std::string(), "Vector of opposite side tags to generate, in the format \033[3m outputTreeName decayDescriptor \033[0m.").getVector();
  std::string BDataFile       = NamedParameter<std::string>("BDataFile", "", "Data file for BDecays");
  std::string BIntFile       = NamedParameter<std::string>("BIntFile", "", "Integration file for BDecays");
  std::string BlogFile      = NamedParameter<std::string> ("BLogFile", "BFitter.log", "Log File for combined B fit");
  auto scanName = NamedParameter<std::string>("scanName", "pCoherentSum::gamma");
  auto scanOutput = NamedParameter<std::string>("scanOutput","scan.txt");
  bool doBScan = NamedParameter<bool>("doBScan", false);
  bool doBFit = NamedParameter<bool>("doBFit", false);

 std::vector<EventList> SigData_B;

 std::vector<EventList> SigInt_B;

 std::vector<EventType> SigType_B;
 
 std::vector<std::string> sfList_B;
 std::vector<int> gammaSigns;
  for (auto& BTag : Btags){
  MinuitParameterSet * MPS_B = new MinuitParameterSet();
  MPS_B->loadFromStream();

    std::stringstream B_log;
    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag, ' ')[1];
    int B_Conj = std::stoi(split(BTag, ' ')[2]);
    int gammaSign = std::stoi(split(BTag,' ')[3]);
    B_log<<B_Name<<"_fit.log";
    std::stringstream B_fit;
    B_fit<<B_Name<<"_plots.root";
    auto B_plotName = B_fit.str();
    auto B_logName = B_log.str();

    std::stringstream BDataNameSS;
    std::stringstream BIntNameSS;

    INFO("BDataFile = "<<BDataFile);
    INFO("BIntFile = "<<BIntFile);

    BDataNameSS<<BDataFile<<":"<<B_Name;
    BIntNameSS<<BIntFile<<":"<<B_Name;
    std::string BDataName = BDataNameSS.str();
    std::string BIntName = BIntNameSS.str();

    EventType BEventType = EventType(pNames);

    if (B_Conj == 1){
      BEventType = BEventType.conj(true);
    }

    EventList sigevents_B = EventList(BDataName, BEventType);
    EventList sigMCevents_B = EventList(BIntName, BEventType);
        
    SigData_B.push_back(sigevents_B);
    SigInt_B.push_back(sigMCevents_B);
    SigType_B.push_back(BEventType);
    sfList_B.push_back(B_Pref);
    gammaSigns.push_back(gammaSign);

    

    auto cs_B = pCoherentSum(BEventType,  *MPS_B, B_Pref , gammaSign);
    cs_B.setEvents(sigevents_B);
    cs_B.setMC(sigMCevents_B);
    cs_B.prepare();
    //auto LL_B2 = make_likelihood( events_B["signal"], events_B["tag"], false, cs_B);
    auto LL_B = make_likelihood( sigevents_B, cs_B);
    Minimiser mini_B = Minimiser(LL_B, MPS_B);

    if (doBScan){
    std::stringstream ss_B;
    ss_B<<"Scan_"<<B_Name<<".txt";
    auto scanOutput_B =ss_B.str();

    std::ofstream scanfile;
    scanfile.open(scanOutput_B, std::ios_base::app);
    auto param = (*MPS_B)[scanName];
    double minimum=param->minInit();
    double maximum=param->maxInit();
    double val = minimum;
    double stepSize = param->stepInit();

    param->setCurrentFitVal(val);   


    while (val < maximum){

    

     param->setCurrentFitVal(val);    


      

      INFO("At "<<val<<" FCN = "<<mini_B.FCN());      
      scanfile<<val<<"\t"<<mini_B.FCN()<<"\n";      

      INFO("Norm = "<<cs_B.norm());
      val += stepSize;
      
      }
      
    scanfile.close();

    }    

if (doBFit){
   FitResult * fr_B = Fit(LL_B, cs_B, sigevents_B, sigMCevents_B, *MPS_B, B_logName, inits, maxAttempts, repeatFits);

      doPlots(B_plotName, cs_B, sigevents_B, sigMCevents_B,nBins);
 //   doPlots(B_plotName, cs_B, sigevents_B, tagevents_B, sigMCevents_B, tagMCevents_B, nBins);
  //  fr_B->writeToFile(B_logName);
}




  delete MPS_B;
  }      
  MinuitParameterSet * MPS_B = new MinuitParameterSet();
  MPS_B->loadFromStream();


  INFO("Doing combined fit for B");
  CombLL combLL_B = CombLL(SigData_B, SigInt_B, SigType_B, *MPS_B, sfList_B, gammaSigns);

  Minimiser combMini_B = Minimiser(combLL_B, MPS_B);
    combMini_B.gradientTest();


     if (doBScan){ 
   

    std::stringstream ss_comb;
    ss_comb<<"Scan_"<<"Comb"<<".txt";
    auto scanOutput_comb =ss_comb.str();

    std::ofstream scanfile;
    scanfile.open(scanOutput_comb, std::ios_base::app);
    auto param = (*MPS_B)[scanName];
    double minimum=param->minInit();
    double maximum=param->maxInit();
    double val = minimum;
    double stepSize = param->stepInit();

    param->setCurrentFitVal(val);   


    while (val < maximum){

    

      param->setCurrentFitVal(val);    


      INFO("At "<<val<<" FCN = "<<combMini_B.FCN());      
      scanfile<<val<<"\t"<<combMini_B.FCN()<<"\n";      
      val += stepSize;
      
      }
      
    scanfile.close();
  
     }  


    if (doBFit){
   // combMini.prepare();
    INFO("Minimising now for combined");
    int attempt = 1;
      combMini_B.doFit(); 
      if (repeatFits){
      if (combMini_B.status() != 0){
        INFO("Didn't seem to get a minimum (returned "<<combMini_B.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
      while (attempt < maxAttempts && combMini_B.status() != 0){
        INFO("Didn't seem to get a minimum (returned "<<combMini_B.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
        combMini_B.doFit();
        attempt++;
      }
      }
      }

    FitResult * fr_B = new FitResult(combMini_B); 
    fr_B->print();
    fr_B->writeToFile(logFile);
    std::map<std::string, std::vector<double> > fits_B = getParams(*MPS_B);
    std::map<std::string, double> pulls_B = getPulls(fits_B, inits);
    for(map<std::string, double >::iterator it = pulls_B.begin(); it != pulls_B.end(); ++it) {
      INFO("Pull = "<<it->first<<" "<<it->second);
    }
    writePulls(BlogFile, pulls_B);
   fr_B->writeToFile(BlogFile);
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
FitResult* Fit( PDF&& ll, pCorrelatedSum& cs, EventList data1, EventList data2, EventList mc1, EventList mc2, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits, int maxAttempts, bool repeatFits )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  MPS.loadFromStream();

  ll.setEvents( data1, data2 );
 
  Minimiser mini( ll, &MPS );
   
  mini.prepare();
  mini.gradientTest();
 INFO("Fitting now"); 
  
  INFO("Making Fit Result");

  int attempt = 1;
  mini.gradientTest();
  mini.doFit();
  if (repeatFits){
  if (mini.status()!=0){

    INFO("Didn't seem to get a minimum (returned "<<mini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
    while (attempt < maxAttempts && mini.status() != 0){
    INFO("Didn't seem to get a minimum (returned "<<mini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);

    mini.doFit();
    attempt++;
    }
  }
  }
  FitResult* fr = new FitResult(mini);



  auto dataEvents = data1;
  auto mcEvents = mc1;
  auto tagMCEvents = mc2;


  int minEvents=15;
  auto binning = BinDT( dataEvents, MinEvents( minEvents ), Dim( dataEvents.eventType().dof() ) );

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
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );

  fr->print();
  return fr;
}


template <typename PDF>
FitResult* Fit( PDF&& ll, pCoherentSum& cs, EventList data1, EventList mc1, MinuitParameterSet& MPS , std::string logFile,   std::map<std::string, std::vector<double> > inits, int maxAttempts, bool repeatFits )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  MPS.loadFromStream();

  ll.setEvents( data1 );
 
  Minimiser mini( ll, &MPS );
   
  mini.prepare();
  mini.gradientTest();
 INFO("Fitting now"); 
  
  INFO("Making Fit Result");

  int attempt = 1;
  mini.gradientTest();
  mini.doFit();
  if (repeatFits){
  if (mini.status()!=0){

    INFO("Didn't seem to get a minimum (returned "<<mini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);
    while (attempt < maxAttempts && mini.status() != 0){
    INFO("Didn't seem to get a minimum (returned "<<mini.status()<<" , trying "<<attempt<<"/"<<maxAttempts);

    mini.doFit();
    attempt++;
    }
  }
  }
  FitResult* fr = new FitResult(mini);



  auto dataEvents = data1;
  auto mcEvents = mc1;



  int minEvents=15;
  auto binning = BinDT( dataEvents, MinEvents( minEvents ), Dim( dataEvents.eventType().dof() ) );

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

    if ( j % 1000000 == 0 && j != 0 ) INFO( "Binned " << j << " sim. events" );
    double w = cs.prob( evt1 ) * (evt1.weight() / evt1.genPdf());
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
        sgn = reOrIm == "Re" ? p.CP() : 1; 
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

template <typename pdftype>
void doPlots(std::string plotFile, pdftype&& cs,  EventList sigEvents, EventList sigMCEvents, int nBins)
{
    

    auto signalType = sigEvents.eventType();

    INFO( "norm[1] = " << cs.norm() );
    TFile* f = TFile::Open(plotFile.c_str(),"RECREATE");

    auto projections = signalType.defaultProjections(nBins);
    for( auto& projection : projections ){
      auto data_plot = projection(sigEvents);
      auto hist = projection.plot();
      for(unsigned i = 0 ; i != sigMCEvents.size(); ++i)
      {
        hist->Fill( projection( sigMCEvents[i] ), cs.prob( sigMCEvents[i]  ) );
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
          hdalitz->Fill( pos.first, pos.second, cs.prob( sigMCEvents[event] ) );
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

    CorrelatedLL<EventList, pCorrelatedSum&> csLL =  make_likelihood(sigEvents, tagEvents, cs);
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
