#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/CombCorrLL.h"
#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/polyLASSO.h"
#include "AmpGen/ProfileClock.h"
#include <TMath.h>
#include "AmpGen/AddCPConjugate.h"

//#include <Math/IFunction.h>
#include <Math/Functor.h>
#include "AmpGen/PhaseCorrection.h"
#include <typeinfo>

//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;




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


EventList getEvents(std::string type, std::vector<std::string> sigName, std::string tagName, std::string dataFile){
   


    INFO("From "<<dataFile);
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


    
    if (type=="signal"){
     return sigEvents;
    }
    else if (type=="tag"){ 
        return tagEvents;
    }
    else {
        return null;
    }
    


}




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

  auto BTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();

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
//    add_CP_conjugate(MPS);
      AddCPConjugate(MPS);
  }
  std::map<std::string, std::vector<double> > inits = getParams(MPS);
//auto fs = std::vector<std::function<double(void)> > {};

  int NInt = NamedParameter<int>("NInt", 1e7);
  INFO("Using "<<NInt<<" integration events");
  EventType eventType(pNames);
  EventList mc =  Generator<>(eventType, &rndm).generate(NInt);
  INFO("Generated Integration events");

  if (NamedParameter<bool>("testPhaseCorr", true)){
    auto event = mc[0];
    INFO("Making phase correction");
    PhaseCorrection phaseCorr(MPS);
    //Expression corr = phaseCorr.calcCorr(event);
    auto corrF = phaseCorr.makefunc();
    size_t order = NamedParameter<size_t>("pCorrelatedSum::Order");
    auto x = event.s(0,1);
    auto y = event.s(0,2);
    auto z = event.s(1,2);
    auto mp = sqrt(event.s(1,1))/2;
    auto mm = sqrt(event.s(2,2))/2;
    auto mK = sqrt(event.s(0,0))/2;
    auto mD = sqrt(x + y + z - pow(mp,2) - pow(mm,2) - pow(mK,2) ) ;
    real_t xmin = pow(mp + mK, 2);
    real_t xmax = pow(mD - mm, 2);
    real_t x0 = (xmax + xmin)/2;
    real_t ymin = pow(mp + mK, 2);
    real_t ymax = pow(mD - mp, 2);
    real_t y0 = (ymax + ymin)/2;
    real_t X = (2 * x - xmax - xmin)/(xmax - xmin);
    real_t Y = (2 * y - ymax - ymin)/(ymax - ymin);
    real_t zp = (X+Y)/2;
    real_t zm = (X-Y)/2;
    INFO("xmin = "<<xmin);
    INFO("xmax = "<<xmax);
    INFO("ymin = "<<ymin);
    INFO("ymax = "<<ymax);

/*
    ProfileClock pcF;
    pcF.start();
    auto F = corrF(zp,zm,order,MPS);
    pcF.stop();
    INFO("stored f("<<zp<<", "<<zm<<") = "<<F<<" took "<<pcF.t_duration);

    ProfileClock pcNormCorr;
    pcNormCorr.start();
    auto normF = phaseCorr.calcCorr(event);
    pcNormCorr.stop();

    INFO("normal f("<<zp<<", "<<zm<<") = "<<normF<<" took "<<pcNormCorr.t_duration);

    ProfileClock pcFastCorr;
    pcFastCorr.start();
    auto fastF = phaseCorr.fastCorr(event);
    pcFastCorr.stop();
    INFO("explicit f("<<zp<<", "<<zm<<") = "<<fastF<<" = "<<fastF()<<" took "<<pcFastCorr.t_duration);
    */
   
    ProfileClock pcCalcCorr;
    pcCalcCorr.start();
    auto corr = phaseCorr.calcCorrL(event);
    pcCalcCorr.stop();
    INFO("calc corr f("<<zp<<", "<<zm<<" ) = "<<corr<<" took "<<pcCalcCorr.t_duration); 



    return 0;
  }


INFO("Doing loop of Fits");

 std::vector<EventList> SigData;
 std::vector<EventList> TagData;
 std::vector<EventList> SigInt;
 std::vector<EventList> TagInt;
 std::vector<EventType> SigType;
 std::vector<EventType> TagType;
 std::vector<std::string> sumFactors;
 std::vector<bool> constAmps;
 SimFit simfit;

std::vector<CorrelatedLL<EventList, pCorrelatedSum&> > pdfs;


for (int i=0; i < tags.size(); i++){
   std::stringstream tag_log;
    auto tagName = split(tags[i],' ')[0];
    tag_log<<tagName<<"_fit.log";
    std::stringstream tag_fit;
    tag_fit<<tagName<<"_plots.root";
    auto tag_plotName = tag_fit.str();
    auto tag_logName = tag_log.str();

    auto sigevents_tag = getEvents("signal", pNames, tags[i], dataFile);
    INFO("Have signal events for "<<i);

    auto tagevents_tag = getEvents("tag", pNames, tags[i], dataFile);

    INFO("Have tag events for "<<i);

    auto types = makeEventTypes(pNames, tags[i]);

    auto signalType = types["signal"];
    auto tagType = types["tag"];
    if (signalType==tagType){
      constAmps.push_back(false);
    }
    else{
      constAmps.push_back(true);
    }
 




    SigData.emplace_back(sigevents_tag);
    TagData.emplace_back(tagevents_tag);
    SigType.emplace_back(sigevents_tag.eventType());
    TagType.emplace_back(tagevents_tag.eventType());
   
   
} 

 //return 0;   
    


 
  std::vector<EventList> BSigData;
  std::vector<std::string> BsumFactors;
  std::vector<int> BgammaSigns;
  std::vector<int> BuseXYs;
  std::vector<int> BConj;



  for (auto& BTag : BTags){
    INFO("B DecayType = "<<BTag); 
    EventType eventType = EventType(pNames);
    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag,' ')[1];
    int B_Conj = std::stoi(split(BTag,' ')[2]);
    int gammaSign = std::stoi(split(BTag,' ')[3]);
    int useXY = std::stoi(split(BTag,' ')[4]);
    
    INFO("GammaSign = "<<gammaSign);


    std::string DataFile = NamedParameter<std::string>("BDataSample", "");
    std::string IntFile = NamedParameter<std::string>("BIntegrationSample", "");

    std::stringstream DataSS;
    DataSS<<DataFile<<":"<<B_Name;
    std::string DataLoc = DataSS.str();


    std::stringstream IntSS;
    IntSS<<IntFile<<":"<<B_Name;
    std::string IntLoc = IntSS.str();

    EventList Data = EventList(DataLoc, eventType);
//    EventList Int = EventList(IntLoc, eventType);


//    auto sig = pCoherentSum(eventType, MPS ,B_Pref, gammaSign, useXY, B_Conj);
//    sig.setEvents(Data);
//    sig.setMC(Int);


    BSigData.emplace_back(Data);
    BgammaSigns.emplace_back(gammaSign);
    BuseXYs.emplace_back(useXY);
    BConj.emplace_back(B_Conj);


  }
  /*
  CombGamCorrLL(std::vector<EventList> SigData, 
		   std::vector<EventList> TagData, 

		   EventType SigType,
		   std::vector<EventType> TagType,
           std::vector<bool> constAmp,

		   std::vector<EventList> GamData, 
           std::vector<int> gammaSigns,
           std::vector<int> useXYs,

           std::vector<int> BConj,
            MinuitParameterSet mps):
  
  */

  bool useCache = NamedParameter<bool>("useCache", true);
  if (useCache){

      CombGamCorrLL LL(SigData, TagData, SigType[0], TagType, constAmps, BSigData, BgammaSigns, BuseXYs, BConj, MPS);
      real_t LL_corr0 = LL.LLCorr(0);
      INFO("LLCorr0 = "<<LL_corr0);
      ProfileClock LLClock;
      LLClock.start();
      real_t ll = LL.LL();
      LLClock.stop();
      INFO("LL = "<<ll<<" took "<<LLClock.t_duration/1000<<" s");
      simfit.add(LL);

      polyLASSO lasso(simfit, MPS);
      Minimiser mini = Minimiser(lasso, &MPS);
      mini.gradientTest();
      ProfileClock miniClock;
      miniClock.start();
      mini.doFit();
      miniClock.stop();
      INFO("Minimiser took "<<miniClock.t_duration/1000<<" s");
      FitResult * fr = new FitResult(mini);
      fr->writeToFile(logFile);
      return 0;
  }
  else{
      EventList mc =  Generator<>(eventType, &rndm).generate(NInt);
      SimFit sfLL;
      for (int i=0; i < SigData.size(); i++){ 
        pCorrelatedSum pdf(SigType[i], TagType[i], MPS);
        EventList tagMC = Generator<>(TagType[i], &rndm).generate(NInt);
        pdf.setEvents(SigData[i], TagData[i]);
        pdf.setMC(mc, tagMC);
        pdf.prepare();
        INFO("pdf(0,0) = "<<pdf.getVal(SigData[i][0], TagData[i][0]));
        real_t norm = pdf.norm();
        INFO("norm = "<<norm);

        INFO("prob(0,0) = "<<pdf.prob(SigData[i][0], TagData[i][0]));
        real_t nA;
        real_t nB;
        real_t nC;
        real_t nD;
        for (size_t k=0;k<SigData[i].size();k++){
          auto v = pdf.getVals(SigData[i][k], TagData[i][k]);
          auto a = v[0];
          auto b = v[1];
          auto c = v[2];
          auto d = v[3];
          nA += std::norm(a);
          nB += std::norm(b);
          nC += std::norm(c);
          nD += std::norm(d);
          INFO("a = "<<a);
          INFO("|a|^2 = "<<std::norm(a));
          INFO("|a|^2 = "<<std::pow(std::abs(a), 2));

        }
        
  
               
        //CorrelatedLL<EventList, pCorrelatedSum&> LLCorr = make_likelihood(SigData[i], TagData[i], pdf);
        
        //INFO("LLCorr = "<<LLCorr.getVal());
        /*
        sfLL.add(LLCorr);
        INFO("sfLL = "<<sfLL.getVal());
        */
      } 
      return 0;

  }
//return 0;

/*
    INFO("Making Combined Minimiser object");
    std::vector<SumPDF<EventList, pCoherentSum&>> pdfsB;


 CombCorrLL corrLL = CombCorrLL(SigData, TagData, SigInt, TagInt, SigType, TagType, MPS, sumFactors);
    simfit.add(corrLL);
    
  pdfsB.reserve(BSigData.size());
  std::vector<pCoherentSum> fcsB(BSigData.size());
  for (size_t i=0;i<BSigData.size(); i++){
   fcsB[i] = pCoherentSum(eventType, MPS, BsumFactors[i], BgammaSigns[i], BuseXYs[i], BConj[i]);
   pdfsB.emplace_back( make_pdf(fcsB[i])); 
   pdfsB[i].setEvents(BSigData[i]);
//   auto& mc  = SigInt[i];
   for_each(pdfsB[i].pdfs(), [&mc](auto& pdf){pdf.setMC(mc);});
   simfit.add(pdfsB[i]);
}
*/
/*
INFO("Making Combined LL");
 auto combLL =  CombGamCorrLL(
        SigData, 
        TagData, 
        BSigData, 
        SigInt, 
        TagInt, 
        BSigInt, 
        SigType,
        TagType,
        BSigType,
        MPS,
        BsumFactors,
        BgammaSigns,
        BuseXYs,
        BConj);

    INFO("CombCorrLL = "<<combLL.getVal());
//    auto commLL2 = SumLL(_LLs);
    simfit.add(combLL);
    INFO("Making Combined Minimiser object");
//    INFO("totalLL = "<<totalLL.getVal());  


    if (NamedParameter<bool>("testLasso", true)){
      INFO("Making LASSO");
    //LASSO lasso(MPS);

      ProfileClock pc;
      pc.start();


      INFO("LL = "<<simfit.getVal());
      pc.stop();
      INFO("time to calc LL = "<<pc.t_duration);
      polyLASSO lasso(simfit, MPS);
      double penTerm = lasso.penalty();
      double lassoVal = lasso.getVal();

      INFO("penalty = "<<penTerm);
      INFO("LL + pen = "<<lassoVal);
      Minimiser lassoMini = Minimiser(lasso, &MPS);


    lassoMini.gradientTest();
      lassoMini.doFit();
        std::ofstream paramsFile;
  paramsFile.open(NamedParameter<std::string>("ParamOutput", "freeParams.opt"));
  for (auto & p:MPS){
	  if (p->isFree()){
		  paramsFile << p->name() << " Free "<<" " << p->mean() << " " << p->err() << "\n ";
		  std::cout << p->name() << " Free "<<" " << p->mean() << " " << p->err() << "\n ";

	  }
	  else {

		  paramsFile << p->name() << " Fix "<<" " << p->mean() << " " << p->err() << "\n ";
	  }
  }
paramsFile.close();



      FitResult * frLasso = new FitResult(lassoMini);
      frLasso->print();
      frLasso->writeToFile(logFile);

      return 0;
    }
    
    Minimiser combMini = Minimiser(simfit, &MPS);
    combMini.gradientTest();
   // combMini.prepare();
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
 std::string covMatrixFile = NamedParameter<std::string>("CovOutput", "gammaCombCorrCov.root");
  auto  covMatrix = combMini.covMatrix();
  TFile * fCov = new TFile(covMatrixFile.c_str(), "recreate");
  fCov->cd();
  covMatrix.Write("CovMatrix"); 
  fCov->Close();
  
*/
  return 0;
}
