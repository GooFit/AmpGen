#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/AddCPConjugate.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"
#include "AmpGen/QMI.h"
#include "AmpGen/Chi2Estimator.h"


using namespace AmpGen;

int main(int argc , char* argv[] ){
  OptionsParser::setArgs( argc, argv );
  MinuitParameterSet MPS;
  MPS.loadFromStream();



  EventType sigType(NamedParameter<std::string>("EventType", "", "Signal Type to generate"));
  size_t NInt(NamedParameter<size_t>("NInt", 1e7, "Number of events to calculate normalisation - should be large"));
  size_t seed(NamedParameter<size_t>("Seed", 1, "Random seed for generation"));
    TRandom3 rndm(seed);
//  size_t nEvents(NamedParameter<size_t>("nEvents", 1000, "number of events to generate, multiplies by the BR of each tag"));
  size_t plot_nBins(NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms"));
  auto tags = NamedParameter<std::string>("TagTypes").getVector();
  auto btags = NamedParameter<std::string>("BTagTypes").getVector();
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  std::string logFile(NamedParameter<std::string>("LogFile", "Fit.log"));
  std::string plotFile(NamedParameter<std::string>("Plots", "Fit.root"));
  std::string psi3770Log(NamedParameter<std::string>("psi3770Log", ""));
  
  #ifdef _OPENMP
    omp_set_num_threads( nThreads );
    INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
    omp_set_dynamic( 0 );
  #endif



  std::string BESIIIFile(NamedParameter<std::string>("BESIIIDataSample", "besiii.root", "Psi(3770) to DDbar events"));
  std::string LHCbFile(NamedParameter<std::string>("LHCbDataSample", "lhcb.root", "B+- to D h+- events"));
  
  EventList_type mcSig =  Generator<>(sigType, &rndm).generate(NInt);
//  EventList_type mcSig2 =  Generator<>(sigType, &rndm).generate(int(1e5));

  std::map<std::string, CompiledExpression<real_t(const real_t*, const real_t*)> > compiledPoly(QMI::cPhaseCorrection(sigType, MPS));

  MinuitParameterSet MPS_Kspipi;
  MPS_Kspipi.loadFromFile("Kspipi.opt");
  AddCPConjugate(MPS_Kspipi);
  CoherentSum A(sigType, MPS_Kspipi);

  CoherentSum Abar(sigType.conj(true), MPS_Kspipi);

  A.setEvents(mcSig);
  A.setMC(mcSig);
  A.prepare();
  Abar.setMC(mcSig);
  Abar.prepare();
  real_t ANorm(A.norm());
  real_t AbarNorm(Abar.norm());
  std::map<std::string, std::vector<real_t> > AampsMC(QMI::AmpArrays(mcSig, A, Abar));
  std::vector<std::map<std::string, std::vector<real_t> > > BampsMC(tags.size());
  std::vector<std::map<std::string, std::vector<real_t> > > AampsDataBESIII(tags.size());
  std::vector<std::map<std::string, std::vector<real_t> > > BampsDataBESIII(tags.size());
  std::vector<EventList> dataSigBESIII(tags.size());
  std::vector<EventList> dataTagBESIII(tags.size());

  std::vector<real_t> BnormsBESIII(tags.size());
  std::vector<real_t> BbarnormsBESIII(tags.size());
  std::vector<real_t> CTagBESIII(tags.size());
  std::vector<real_t> STagBESIII(tags.size());
  std::vector<std::string> BESIIINames(tags.size());

  std::vector<bool> BESIIISameType(tags.size());



  std::vector<std::string> LHCbNames(btags.size());
  std::vector<EventList> dataLHCb(btags.size());
  std::vector<std::map<std::string, std::vector<real_t> > > AampsDataLHCb(btags.size());
  std::vector<Int_t> gammaSigns(btags.size());

  for(unsigned i=0;i<btags.size();i++){
    auto a = split(btags[i], ' ');
    Int_t gammaSign(std::stoi(a[1]));
    gammaSigns[i] = gammaSign;
    LHCbNames[i] = a[0];
    EventList data(LHCbFile + ":" + a[0], sigType);
    dataLHCb[i] = data;
    complex_t zB = QMI::ckm_zB(MPS, gammaSign);
    INFO("rB = "<<std::abs(zB)<<" thetaB = "<<std::arg(zB));


    //std::map<std::string, std::vector<real_t> > ampsDat(QMI::AmpArrays(data, A, Abar));
    AampsDataLHCb[i] = QMI::AmpArrays(data, A, Abar);
//
//    std::function<real_t()> myLL = [&ANorm, &AbarNorm, &gammaSign, &MPS, &mcSig, &ampsDat, &ampsMC, &compiledPoly, &data](){
//        complex_t zB(QMI::ckm_zB(MPS, gammaSign));
//        real_t C(QMI::cosTerm(ampsMC, MPS, compiledPoly, mcSig));
//        real_t S(QMI::sinTerm(ampsMC, MPS, compiledPoly, mcSig));
//        real_t norm = 0;
//        if (gammaSign>0) norm = AbarNorm + std::norm(zB) * ANorm - 2 * std::real(zB)*C - 2*gammaSign * std::imag(zB) * S;
//        if (gammaSign<0) norm = ANorm + std::norm(zB) * AbarNorm - 2 * std::real(zB)*C - 2*gammaSign * std::imag(zB) * S;
//        real_t LL =0;
//        #pragma omp parallel for reduction (+:LL)
//        for (unsigned i=0;i<data.size();++i){
//          real_t correction = QMI::phaseCorrection(data[i], compiledPoly, MPS);
//          real_t prob = 0;
//          if (gammaSign>0)  prob = std::norm(ampsDat["Abar"][i]) + std::norm(zB * ampsDat["A"][i]) + 2 * std::abs(zB) * ampsDat["A"][i] * ampsDat["Abar"][i] * cos(ampsDat["dd"][i] + correction - gammaSign * std::arg(zB));
//          if (gammaSign<0)  prob = std::norm(ampsDat["A"][i]) + std::norm(zB * ampsDat["Abar"][i]) + 2 * std::abs(zB) * ampsDat["A"][i] * ampsDat["Abar"][i] * cos(ampsDat["dd"][i] + correction - gammaSign * std::arg(zB));
//          LL += log(prob/norm); 
//        }
//      return -2 * LL;
//    };//QMI::LLCKM(data, MPS, A, Abar, gammaSign, compiledPoly, mcSig);

    //INFO("Made LL "<<a[0]<<" = "<<LLLHCb[i]());
  }


  auto LL_LHCb = [&dataLHCb, &mcSig, &gammaSigns, &AampsDataLHCb, &AampsMC, &MPS, &compiledPoly, &ANorm, &AbarNorm](){
    real_t ll=0;
      real_t C = QMI::cosTerm(AampsMC, MPS, compiledPoly, mcSig);
      real_t S = QMI::sinTerm(AampsMC, MPS, compiledPoly, mcSig);
    for (unsigned i =0;i<dataLHCb.size();++i){
      complex_t zB(QMI::ckm_zB(MPS, gammaSigns[i]));
      real_t norm =0 ;
      real_t normNonInt = 0;
      real_t normInt = 0;

      if (gammaSigns[i]>0){
        normNonInt = ANorm * std::norm(zB) + AbarNorm;
      }
      else{
        normNonInt = ANorm + AbarNorm * std::norm(zB);
      }
      normInt = 2 * C * std::real(zB) -2 * gammaSigns[i] * S * std::imag(zB);
      norm = normNonInt + normInt;
 //     INFO("norm = "<<norm);

      real_t ll_tag =0;
      #pragma omp parallel for reduction (+:ll_tag)
      for (unsigned j=0;j<dataLHCb[i].size();++j){
        real_t prob = 0;
        real_t correction = QMI::phaseCorrection(dataLHCb[i][j], compiledPoly, MPS);
        if (gammaSigns[i]>0){
          prob = std::norm(AampsDataLHCb[i]["A"][j] * zB) + std::norm(AampsDataLHCb[i]["Abar"][j]) + 2 * AampsDataLHCb[i]["A"][j] * AampsDataLHCb[i]["Abar"][j] * std::abs(zB) * cos(AampsDataLHCb[i]["dd"][j] + correction - gammaSigns[i] * std::arg(zB));
        }
        else{
          prob = std::norm(AampsDataLHCb[i]["A"][j]) + std::norm(AampsDataLHCb[i]["Abar"][j]  * zB) + 2 * AampsDataLHCb[i]["A"][j] * AampsDataLHCb[i]["Abar"][j] * std::abs(zB) * cos(AampsDataLHCb[i]["dd"][j] + correction - gammaSigns[i] * std::arg(zB));

        }
//        INFO(j<<" "<<" prob = "<<prob);
        ll_tag += log(prob/norm);
      }
      ll += ll_tag;
    } 
    return -2 * ll;
  };

  ProfileClock  tLHCb;

  tLHCb.start();
  real_t LL_LHCb_0(LL_LHCb());
  tLHCb.stop();
  INFO("Start with "<<LL_LHCb_0<<" took "<<tLHCb.t_duration<<"ms to calculate LL_LHCb");
  
//  auto LL_BESIII_Total = [&LLBESIII](){
//    real_t ll = 0;
//    for (unsigned i=0;i<LLBESIII.size();++i){
//      INFO(LLBESIII[i]());
//      ll += LLBESIII[i]();
//    }
//    return ll;
//  };
//
//  auto LL_LHCb_Total = [&LLLHCb](){
//    real_t ll = 0;
//    for (unsigned i=0;i<LLLHCb.size();++i){
//
//      INFO(LLLHCb[i]());
//      ll += LLLHCb[i]();
//    }
//    return ll;
//  };
//  INFO("LL = "<<LL_BESIII_Total());
//  INFO("LL = "<<LL_LHCb_Total());
 // INFO("LL = "<<LL_BESIII_LHCb_Total());


  FitResult* fr;
  TFitResult * root_fr;
  TFile * tOutFile = TFile::Open(plotFile.c_str(), "RECREATE");
  tOutFile->cd();
 bool doFit = NamedParameter<bool>("doFit", true);   
 bool writeTuple = NamedParameter<bool>("writeTuple", true);   

 if (psi3770Log != ""){
   INFO("Getting mu/sigma for parameters");
   auto meanAndErr_psi3770 = QMI::fitValAndErr(psi3770Log);
   INFO("Getting inverse Covariance matrix");
   auto invCovMatrix_psi3770 = QMI::invCovarianceMatrix(psi3770Log);
   INFO("Done invConv");
   auto constraint = [&meanAndErr_psi3770, &invCovMatrix_psi3770, &MPS](){
     return QMI::myGaussConstraint(meanAndErr_psi3770, invCovMatrix_psi3770, MPS);
   };
   INFO("Built Constraint as function");

   auto LL_LHCb_Constrained = [&LL_LHCb, &constraint](){
     return LL_LHCb() + constraint();
   };

//   mini = Minimiser(LL_LHCb_Constrained, &MPS); 
    Minimiser mini(LL_LHCb_Constrained, &MPS);
    
    if (doFit){ 
        
mini.gradientTest();
mini.doFit();
    fr = new FitResult(mini);
    
    root_fr = new TFitResult(mini.fitResult());
   
    }


    
//    fr.print();
 //   fr.writeToFile(logFile);
 //   fr(_fr);

 }
  else{
    Minimiser mini(LL_LHCb, &MPS);
    if (doFit){
mini.gradientTest();
        mini.doFit();

  



//    fr.print();
//    fr.writeToFile(logFile);
//    fr(_fr);
  }
fr = new FitResult (mini);
    root_fr = new TFitResult(mini.fitResult());
  }
  


    root_fr->SetName("FitResult");
    root_fr->Write();
    
  std::vector<real_t> chi2(dataLHCb.size());
  std::vector<real_t> nBins(dataLHCb.size());

  auto my_dd = [&A, &Abar](Event& evt){
    return QMI::dd(evt, A, Abar);
  };
  auto my_corr = [&compiledPoly, &MPS](Event& evt){
    return QMI::phaseCorrection(evt, compiledPoly,  MPS);
  };

  QMI::writeValues(mcSig, my_dd, "dd");
  QMI::writeValues(mcSig, my_corr, "corr");
  QMI::writeDalitz(mcSig);

  size_t i=0;
  real_t totalChi2 = 0;
  size_t totalnBins = 0;
  real_t C = QMI::cosTerm(AampsMC, MPS, compiledPoly, mcSig);
  real_t S = QMI::sinTerm(AampsMC, MPS, compiledPoly, mcSig);

  for (auto tag : btags){
      INFO("tag = "<<tag);
      auto a = split(tag, ' ');
//      EventType tagType(Particle(a[1], {}, false).eventType());
//      EventList_type mcTag =  Generator<>(tagType, &rndm).generate(NInt);
      INFO("a[0] = "<<a[0]);
//      size_t nEvents_tag(std::stod(a[2]) * nEvents);

//      bool sameType = sigType == tagType;
      
//      real_t C = QMI::cosTermNoCorrection(AampsMC, mcSig);
//      real_t S = QMI::sinTermNoCorrection(AampsMC,  mcSig);
     real_t _chi2=0;
      size_t _nBins=0;
     const std::function<double(const Event&)> psi = [&A, &Abar, &MPS, &compiledPoly, &gammaSigns, &i](const Event& evt){
          return QMI::probCKM_unnorm(evt, A, Abar, MPS, compiledPoly, gammaSigns[i]);///norm;
    };
    INFO("Doing Chi2");
    Chi2Estimator chi2Est(dataLHCb[i], mcSig, psi, MinEvents(NamedParameter<size_t>("MinEvents", 1)), Dim(dataLHCb[i].eventType().dof()) );
    _chi2 = chi2Est.chi2();
    _nBins = chi2Est.nBins();
      int nFree = 0;
      for (auto& p : MPS){
          if (p->isFree()) nFree += 1;
      }
   //   auto chi2_par = MPS.add(a[0] + "_chi2", Flag::Fix, _chi2, 0);
   //   auto ndf_par = MPS.add(a[0] + "_dof", Flag::Fix, _nBins - nFree - 1,0);
          
   
             TVectorD tvchi2(1);
          TVectorD tvndf(1);
          tvchi2[0] = _chi2;
          tvndf[0] = _nBins - nFree - 1;
          tvchi2.Write((a[0] + "_chi2").c_str());
          tvndf.Write((a[0] + "_ndf").c_str());




    chi2[i] = _chi2;
    nBins[i] = _nBins;
    totalChi2 += _chi2;
    totalnBins += _nBins;

    INFO("Doing Data Proj");
    auto data_projection_s01 = dataLHCb[i].eventType().projection(plot_nBins, {0, 1});
    auto data_projection_s02 = dataLHCb[i].eventType().projection(plot_nBins, {0, 2});
    auto data_projection_s01_vs_s02 = dataLHCb[i].makeProjection(Projection2D(data_projection_s01, data_projection_s02), PlotOptions::Prefix(("Data_s01_vs_s02_" + a[0]).c_str()));
    data_projection_s01_vs_s02->Write();
    TH1D* data_hist_s01 = dataLHCb[i].makeProjection(data_projection_s01, PlotOptions::Prefix(("Data_s01_" + a[0]).c_str()));
    TH1D* data_hist_s02 = dataLHCb[i].makeProjection(data_projection_s02, PlotOptions::Prefix(("Data_s02_" + a[0]).c_str()));
    data_hist_s01->Write();
    data_hist_s02->Write();
    INFO("Doing MC Proj");
    ArgumentPack s01_args_fit(PlotOptions::Prefix("Fit_s01_" + a[0]),
                              PlotOptions::Norm(dataLHCb[i].size()),
                              WeightFunction(psi));
    ArgumentPack s02_args_fit(PlotOptions::Prefix("Fit_s02_" + a[0]),
                              PlotOptions::Norm(dataLHCb[i].size()),
                              WeightFunction(psi));
    ArgumentPack s01_vs_s02_args_fit(PlotOptions::Prefix("Fit"  + a[0]),
          QMI::QMIPlotOptions::name("s01_vs_s02_" + a[0]),
          QMI::QMIPlotOptions::posXFunction([](Event evt){return evt.s(0, 1);} ),
          QMI::QMIPlotOptions::posYFunction([](Event evt){return evt.s(0, 2);} ),
          PlotOptions::Bins(plot_nBins),
          PlotOptions::Norm(dataLHCb[i].size())
           );



    auto mc_plot_s01 = data_projection_s01(mcSig, s01_args_fit);
    auto mc_plot_s02 = data_projection_s02(mcSig, s02_args_fit);
    auto mc_plot_s01_vs_s02 = QMI::proj_2D_BDK(mcSig, psi, s01_vs_s02_args_fit);
    mc_plot_s01->Write();
    mc_plot_s02->Write();
    mc_plot_s01_vs_s02->Write();
    i++;
  }
  fr->addChi2(totalChi2, totalnBins);
  fr->print();
  fr->writeToFile(logFile);
  
  tOutFile->Close();


   
//

  return 0;
}

//template <typename=pdf_t> void doFit()
