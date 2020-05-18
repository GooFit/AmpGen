#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/CombCorrLL.h"
//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

std::vector<std::string> makeBranches(EventType Type, std::string prefix);
template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data1, EventList& data2, EventList& mc1, EventList& mc2, MinuitParameterSet& MPS );

MinuitParameterSet *  copyMPS(MinuitParameterSet& mps);

std::map<std::string, double>  getPulls(std::map<std::string, std::vector<double> > fit, std::map<std::string, std::vector<double> > init);

void writePulls(std::string fileName, std::map<std::string, double>);
std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps);


EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile, std::string intFile);




std::map<std::string, EventType> makeEventTypes(std::vector<std::string> sigName, std::string tagName);



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
    int nFits = NamedParameter<int>("nFits", 4, "number of repeats of mini.doFits() for debug purposes!");
    bool doProjections = NamedParameter<bool>("doProjections", true);
    bool doPCorrSum = NamedParameter<bool>("doPCorrSum", false);
    auto scanName = NamedParameter<std::string>("scanName", "pCorrelatedSum::C00");
    auto scanOutput = NamedParameter<std::string>("scanOutput","scan.txt");

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "QcFitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");

  bool makeCPConj      = NamedParameter<bool>("makeCPConj", false, "Make CP Conjugates");
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
auto findParam = MPS.find(scanName);
 if (findParam==nullptr) {
    INFO("AmpGen didn't find "<<scanName<<" did you put it in your options file?");
    return 0;
    
 }
 std::vector<EventList> SigData;
 std::vector<EventList> TagData;
 std::vector<EventList> SigInt;
 std::vector<EventList> TagInt;
 std::vector<EventType> SigType;
 std::vector<EventType> TagType;
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
   
   SigData.push_back(sigEvents);
   TagData.push_back(tagEvents);
   SigInt.push_back(sigMCEvents);
   TagInt.push_back(tagMCEvents);
   SigType.push_back(signalType);
   TagType.push_back(tagType);
 
    
    pCorrelatedSum cs(signalType, tagType, MPS);
    cs.setEvents(sigEvents, tagEvents);
    cs.setMC(sigMCEvents, tagMCEvents);
     cs.prepare();  
    auto csLL = make_likelihood(sigEvents, tagEvents, cs);
    csLL.setEvents(sigEvents, tagEvents);
    

 //   INFO( "uncorrected norm = " << cs.norm() );
//    INFO( "corrected norm = " << cs.slowNorm() );
    Minimiser mini( csLL, &MPS );
   
    
 

    std::stringstream ss_tag;
    ss_tag<<"Scan_"<<tokens[0]<<".txt";
    auto scanOutput_tag =ss_tag.str();

    std::ofstream scanfile;
    scanfile.open(scanOutput_tag, std::ios_base::app);
    auto param = MPS[scanName];
    double minimum=param->minInit();
    double maximum=param->maxInit();
    double val = minimum;
    double stepSize = param->stepInit();

    MPS[scanName]->setCurrentFitVal(val);   


    while (val < maximum){

    

     MPS[scanName]->setCurrentFitVal(val);    


      

      INFO("At "<<val<<" FCN = "<<mini.FCN());      
      scanfile<<val<<"\t"<<mini.FCN()<<"\n";      

      INFO("Norm = "<<cs.norm());
      val += stepSize;
      
      }
      
    scanfile.close();
 } 

  CombCorrLL combLL = CombCorrLL(SigData, TagData, SigInt, TagInt, SigType, TagType, MPS);
  INFO("Making Combined Minimiser object");
  bool doComb = true;
  if (doComb){

    
    Minimiser mini_Comb = Minimiser(combLL, &MPS);

    std::stringstream ss_comb;
    ss_comb<<"Scan_"<<"Comb"<<".txt";
    auto scanOutput_comb =ss_comb.str();

    std::ofstream scanfile;
    scanfile.open(scanOutput_comb, std::ios_base::app);
    auto param = MPS[scanName];
    double minimum=param->minInit();
    double maximum=param->maxInit();
    double val = minimum;
    double stepSize = param->stepInit();

    MPS[scanName]->setCurrentFitVal(val);   


    while (val < maximum){

    

      MPS[scanName]->setCurrentFitVal(val);    


      INFO("At "<<val<<" FCN = "<<mini_Comb.FCN());      
      scanfile<<val<<"\t"<<mini_Comb.FCN()<<"\n";      
      val += stepSize;
      
      }
      
    scanfile.close();
  
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
  mini.gradientTest();
  
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


