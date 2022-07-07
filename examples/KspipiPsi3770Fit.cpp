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
#include "AmpGen/QMI.h"
#include "TFitResult.h"

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



  for(unsigned i=0;i<tags.size();++i ){
    auto a = split(tags[i], ' ');
    EventType tagType(Particle(a[1], {}, false).eventType());
    bool sameType = sigType == tagType;
    BESIIISameType[i] = sameType;
    EventList dataSig(BESIIIFile + ":Signal_" + a[0], sigType);
    EventList dataTag(BESIIIFile + ":Tag_" + a[0], tagType);
    //std::map<std::string, std::vector<real_t> > ampsSigDat(QMI::AmpArrays(dataSig, A, Abar));
    AampsDataBESIII[i] = QMI::AmpArrays(dataSig, A, Abar);
    dataSigBESIII[i] = dataSig;
    dataTagBESIII[i] = dataTag;
    BESIIINames[i] = a[0];

    if (sameType){


      //std::map<std::string, std::vector<real_t> > ampsTagDat(QMI::AmpArrays(dataTag, A, Abar));
      BampsDataBESIII[i] = QMI::AmpArrays(dataTag, A, Abar);
      BnormsBESIII[i] = 0;
      BbarnormsBESIII[i] = 0;
      CTagBESIII[i] = 0;
      STagBESIII[i] = 0;
      /*
      std::function<real_t()> myLL = [&ANorm, &AbarNorm, &MPS, &mcSig,  &ampsMC, &compiledPoly, &ampsSigDat, &ampsTagDat, &dataSig, &dataTag](){
        real_t C(QMI::cosTerm(ampsMC, MPS, compiledPoly, mcSig));
        real_t S(QMI::sinTerm(ampsMC, MPS, compiledPoly, mcSig));
        real_t norm = 2 * ANorm * AbarNorm - 2 * C* C - 2  * S *S;

        real_t LL=0;
        #pragma omp parallel for reduction (+:LL)
        for (unsigned i=0;i<dataSig.size();++i){
          real_t correction1 = QMI::phaseCorrection(dataSig[i], compiledPoly, MPS);
          real_t correction2 = QMI::phaseCorrection(dataTag[i], compiledPoly, MPS);
          real_t prob = std::norm(ampsSigDat["A"][i] * ampsTagDat["Abar"][i]) + std::norm(ampsSigDat["Abar"][i] * ampsTagDat["A"][i]) - 2 * ampsSigDat["A"][i] * ampsSigDat["Abar"][i] * ampsTagDat["A"][i] * ampsTagDat["Abar"][i] * cos(ampsSigDat["dd"][i] + correction1 - ampsTagDat["dd"][i] - correction2);
          LL += log(prob/norm);
        }
 
        return -2 * LL;
      };//QMI::LLPsi3770(dataSig, dataTag, MPS, A, Abar, A, Abar, sameType, compiledPoly, mcSig, mcSig);
      */
//      INFO("Made LL "<<a[0]<<" = "<<LLBESIII[i]());
    }
    else{
      MinuitParameterSet MPS_Tag;
      MPS_Tag.loadFromFile(a[0] + ".opt");
      CoherentSum B(tagType, MPS_Tag);
      CoherentSum Bbar(tagType.conj(true), MPS_Tag);

      EventList_type mcTag = Generator<>(tagType, &rndm).generate(NInt);
      B.setMC(mcTag); 
      B.prepare();
      Bbar.setMC(mcTag);
      Bbar.prepare();
//      real_t BNorm(B.norm());
//      real_t BbarNorm(Bbar.norm());
      BnormsBESIII[i] = B.norm();
      BbarnormsBESIII[i] = Bbar.norm();
      //std::map<std::string, std::vector<real_t> > ampsTagMC(QMI::AmpArrays(mcTag, B, Bbar));
      BampsMC[i] = QMI::AmpArrays(mcTag, B, Bbar);
      //std::map<std::string, std::vector<real_t> > ampsTagDat(QMI::AmpArrays(dataTag, B, Bbar));
      BampsDataBESIII[i] = QMI::AmpArrays(dataTag, B, Bbar);
      real_t CTag(QMI::cosTermNoCorrection(BampsMC[i], mcTag));
      real_t STag(QMI::cosTermNoCorrection(BampsMC[i], mcTag));
      CTagBESIII[i] = CTag;
      STagBESIII[i] = STag;

//      INFO("Making LL "<<a[0]);
//      std::function<real_t()> myLL = [&ANorm, &AbarNorm, &BNorm, &BbarNorm, &MPS, &mcSig, &mcTag, &ampsMC, &ampsTagMC, &compiledPoly, &ampsSigDat, &ampsTagDat, &CTag, &STag, &dataSig, &dataTag](){
//        real_t C(QMI::cosTerm(ampsMC, MPS, compiledPoly, mcSig));
//        real_t S(QMI::sinTerm(ampsMC, MPS, compiledPoly, mcSig));
//        real_t norm = ANorm * BbarNorm + AbarNorm * BNorm - 2 * C*CTag - 2 * S * STag;
//        real_t LL=0;
//        #pragma omp parallel for reduction (+:LL)
//        for (unsigned i=0;i<dataSig.size();++i){
//          real_t correction = QMI::phaseCorrection(dataSig[i], compiledPoly, MPS);
//          real_t prob = std::norm(ampsSigDat["A"][i] * ampsTagDat["Abar"][i]) + std::norm(ampsSigDat["Abar"][i] * ampsTagDat["A"][i]) - 2 * ampsSigDat["A"][i] * ampsSigDat["Abar"][i] * ampsTagDat["A"][i] * ampsTagDat["Abar"][i] * cos(ampsSigDat["dd"][i] + correction - ampsTagDat["dd"][i]);
//          LL += log(prob/norm);
//        }
//        return -2* LL;
//        };
      //QMI::LLPsi3770(dataSig, dataTag, MPS, A, Abar, B, Bbar, sameType, compiledPoly, mcSig, mcTag);
//      INFO("Made LL "<<a[0]<<" = "<<LLBESIII[i]());
    }
  }
  

  auto LL_BESIII = [&dataSigBESIII, &dataTagBESIII, &mcSig, &BESIIISameType, &AampsDataBESIII, &BampsDataBESIII, &AampsMC, &BampsMC, &CTagBESIII, &STagBESIII, &MPS, &compiledPoly, &ANorm, &AbarNorm, &BnormsBESIII, &BbarnormsBESIII](){
    real_t ll=0;
      real_t C = QMI::cosTerm(AampsMC, MPS, compiledPoly, mcSig);
      real_t S = QMI::sinTerm(AampsMC, MPS, compiledPoly, mcSig);
    for (unsigned i =0;i<dataSigBESIII.size();++i){
      real_t norm =0 ;
      real_t normNonInt = 0;
      real_t normInt = 0;

      if (BESIIISameType[i]){
        normNonInt =2 * ANorm * AbarNorm;
        normInt = - 2*C *C - 2 * S * S;
      }
      else{
        normNonInt = ANorm * BbarnormsBESIII[i] + AbarNorm * BnormsBESIII[i];
        normInt = -2 * C * CTagBESIII[i] - 2 * S * STagBESIII[i];
      }
      norm = normNonInt + normInt;
      real_t ll_tag =0;
      #pragma omp parallel for reduction (+:ll_tag)
      for (unsigned j=0;j<dataSigBESIII[i].size();++j){
        real_t prob = 0;
        real_t correction1 = QMI::phaseCorrection(dataSigBESIII[i][j], compiledPoly, MPS);
        if(BESIIISameType[i]){
          real_t correction2 = QMI::phaseCorrection(dataTagBESIII[i][j], compiledPoly, MPS);
          prob = std::norm(AampsDataBESIII[i]["A"][j] * BampsDataBESIII[i]["Abar"][j]) + std::norm(AampsDataBESIII[i]["Abar"][j] * BampsDataBESIII[i]["A"][j]) - 2 * AampsDataBESIII[i]["A"][j] * AampsDataBESIII[i]["Abar"][j] * BampsDataBESIII[i]["A"][j] * BampsDataBESIII[i]["Abar"][j] * cos(AampsDataBESIII[i]["dd"][j] + correction1 - BampsDataBESIII[i]["dd"][j] - correction2);
        }
        else{
          prob = std::norm(AampsDataBESIII[i]["A"][j] * BampsDataBESIII[i]["Abar"][j]) + std::norm(AampsDataBESIII[i]["Abar"][j] * BampsDataBESIII[i]["A"][j]) - 2 * AampsDataBESIII[i]["A"][j] * AampsDataBESIII[i]["Abar"][j] * BampsDataBESIII[i]["A"][j] * BampsDataBESIII[i]["Abar"][j] * cos(AampsDataBESIII[i]["dd"][j] + correction1 - BampsDataBESIII[i]["dd"][j]);
        }
        ll_tag += log(prob/norm);
      }
      ll += ll_tag;
    }
    return -2 * ll;
  };

  ProfileClock tBESIII;
  tBESIII.start();
  real_t LL_BESIII_0(LL_BESIII()); 
  tBESIII.stop();
  INFO("Start with "<<LL_BESIII_0<<" took "<<tBESIII.t_duration<<"ms to calculate LL_BESIII");
 
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


     
  Minimiser mini(LL_BESIII, &MPS);
  mini.gradientTest();
  mini.doFit();
  FitResult fr(mini);

  std::vector<real_t> chi2(dataSigBESIII.size());
  std::vector<real_t> nBins(dataSigBESIII.size());
  real_t totalChi2 = 0;
  real_t totalnBins = 0;
  TFile * tOutFile = TFile::Open(plotFile.c_str(), "RECREATE");
  tOutFile->cd();
  size_t i=0;
  for (auto tag : tags){
      INFO("tag = "<<tag);
      auto a = split(tag, ' ');
      EventType tagType(Particle(a[1], {}, false).eventType());
      EventList_type mcTag =  Generator<>(tagType, &rndm).generate(NInt);
//      size_t nEvents_tag(std::stod(a[2]) * nEvents);

      bool sameType = sigType == tagType;
      real_t C = QMI::cosTerm(AampsMC, MPS, compiledPoly, mcSig);
      real_t S = QMI::sinTerm(AampsMC, MPS, compiledPoly, mcSig);
//      real_t C = QMI::cosTermNoCorrection(AampsMC, mcSig);
//      real_t S = QMI::sinTermNoCorrection(AampsMC,  mcSig);
 
      real_t _chi2=0;
      size_t _nBins=0;

      if (!sameType){
        MinuitParameterSet MPS_Tag;
        INFO("Looking for "<<a[0] <<".opt");
        MPS_Tag.loadFromFile(a[0] + ".opt");
      
        CoherentSum B(tagType, MPS_Tag);

        CoherentSum Bbar(tagType.conj(true), MPS_Tag);

        B.setMC(mcTag);
        B.prepare();
        Bbar.setMC(mcTag);
        Bbar.prepare();

        std::map<std::string, std::vector<real_t> > BAmps (QMI::AmpArrays(mcTag, B, Bbar));
        real_t CB = QMI::cosTermNoCorrection(BAmps, mcTag);
        real_t SB = QMI::sinTermNoCorrection(BAmps, mcTag);

        
        //INFO("nB = "<<B.norm());
//        real_t norm(QMI::corr_norm(A, Abar, B, Bbar, MPS, mcSig, compiledPC, sameType));

        real_t norm = A.norm() * Bbar.norm() + Abar.norm() * B.norm() - 2 * C * CB - 2 * S * SB;
        //norm = norm * mcSig.size();

        INFO(a[0]<<" Norm term = "<<norm);
        //real_t probCorr_unnorm(Event& evt1, Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, ce& x, ce& y, bool sameTag)
        auto psi = [&A, &Abar, &B, &Bbar, &MPS, &compiledPoly, &sameType, &norm](Event evt1, Event evt2){
          return QMI::probCorr_unnorm(evt1, evt2, A, Abar, B, Bbar, MPS, compiledPoly, sameType);///norm;
//            return std::norm(A.getValNoCache(evt1) * Bbar.getValNoCache(evt2) - Abar.getValNoCache(evt1) * B.getValNoCache(evt2));
//            return std::norm(A.getValNoCache(evt1));// * Bbar.getValNoCache(evt2) - Abar.getValNoCache(evt1) * B.getValNoCache(evt2));
          };
          QMI::do_chi2_corr(dataSigBESIII[i], dataTagBESIII[i], mcSig, mcTag, psi, MPS, _chi2, _nBins);
          int nFree = 0;
          for (auto& p : MPS){
              if (p->isFree()) nFree += 1;
          }
          //auto chi2_par = MPS.add(a[0] + "_chi2", Flag::Fix, _chi2, 0);
          //auto ndf_par = MPS.add(a[0] + "_dof", Flag::Fix, _nBins - nFree - 1,0);
          
          
         


          ArgumentPack s01_args_fit(PlotOptions::Prefix("Fit" + a[0]),
          QMI::QMIPlotOptions::name("s01_" + a[0]),
          QMI::QMIPlotOptions::posXFunction([](Event evt){return evt.s(0, 1);} ),
          PlotOptions::Bins(plot_nBins),
          PlotOptions::Norm(dataSigBESIII[i].size())
           );
          ArgumentPack s01_vs_s02_args_fit(PlotOptions::Prefix("Fit"  + a[0]),
          QMI::QMIPlotOptions::name("s01_vs_s02_" + a[0]),
          QMI::QMIPlotOptions::posXFunction([](Event evt){return evt.s(0, 1);} ),
          QMI::QMIPlotOptions::posYFunction([](Event evt){return evt.s(0, 2);} ),
          PlotOptions::Bins(plot_nBins),
          PlotOptions::Norm(dataSigBESIII[i].size())
           );

         INFO("Start data Projections");
         auto data_projection_s01 = dataSigBESIII[i].eventType().projection(plot_nBins, {0, 1});
         auto data_projection_s02 = dataSigBESIII[i].eventType().projection(plot_nBins, {0, 2});
         auto data_projection_s01_vs_s02 = dataSigBESIII[i].makeProjection(Projection2D(data_projection_s01, data_projection_s02), PlotOptions::Prefix(("Data_s01_vs_s02_" + a[0]).c_str()));
         INFO("Done data Projections");
          auto mc_projection_s01 = QMI::proj_1D_psi3770(mcSig, mcTag, psi, s01_args_fit);
          auto mc_projection_s01_vs_s02 = QMI::proj_2D_psi3770(mcSig, mcTag, psi, s01_vs_s02_args_fit);
          
          //auto histProj = QMI::proj_1D_psi3770(dataSigBESIII[i], dataTagBESIII[i], psi, PlotOptions::Prefix("Data"));
          INFO("_chi2 = "<<_chi2);
          INFO("_nBins = "<<_nBins);
          totalChi2 += _chi2;
          totalnBins += _nBins;
          TH1D* data_hist_s01 = dataSigBESIII[i].makeProjection(data_projection_s01, PlotOptions::Prefix(("Data_s01_" + a[0]).c_str()));
          TH1D* data_hist_s02 = dataSigBESIII[i].makeProjection(data_projection_s02, PlotOptions::Prefix(("Data_s02_" + a[0]).c_str()));

          data_hist_s01->Write();
          data_hist_s02->Write();
 
          data_projection_s01_vs_s02->Write();
          mc_projection_s01->Write();
          mc_projection_s01_vs_s02->Write();


        }
        else{
//          INFO(a[0]<<" Norm term = "<<norm);

          real_t norm = A.norm() * Abar.norm() + Abar.norm() * A.norm() - 2 * C * C - 2 * S * S;
          auto psi = [&A, &Abar, &MPS, &compiledPoly, &sameType, &norm](Event evt1, Event evt2){
            return QMI::probCorr_unnorm(evt1, evt2, A, Abar, A, Abar, MPS, compiledPoly, sameType)/norm;
            //return std::norm(A.getValNoCache(evt1) * Abar.getValNoCache(evt2) - Abar.getValNoCache(evt1) * A.getValNoCache(evt2));
//            return std::norm(A.getValNoCache(evt1));// * Bbar.getValNoCache(evt2) - Abar.getValNoCache(evt1) * B.getValNoCache(evt2));
            
          };
          QMI::do_chi2_corr(dataSigBESIII[i], dataTagBESIII[i], mcSig, mcTag, psi, MPS, _chi2, _nBins);
          ArgumentPack s01_args_fit(PlotOptions::Prefix("Fit" + a[0]),
          QMI::QMIPlotOptions::name("s01_" + a[0]),
          QMI::QMIPlotOptions::posXFunction([](Event evt){return evt.s(0, 1);} ),
          PlotOptions::Bins(plot_nBins),
          PlotOptions::Norm(dataSigBESIII[i].size())
           );
          ArgumentPack s01_vs_s02_args_fit(PlotOptions::Prefix("Fit" + a[0]),
          QMI::QMIPlotOptions::name("s01_vs_s02_" + a[0]),
          QMI::QMIPlotOptions::posXFunction([](Event evt){return evt.s(0, 1);} ),
          QMI::QMIPlotOptions::posYFunction([](Event evt){return evt.s(0, 2);} ),
          PlotOptions::Bins(plot_nBins),
          PlotOptions::Norm(dataSigBESIII[i].size())
           );

         INFO("Start data Projections");
         auto data_projection_s01 = dataSigBESIII[i].eventType().projection(plot_nBins, {0, 1});
         auto data_projection_s02 = dataSigBESIII[i].eventType().projection(plot_nBins, {0, 2});
         auto data_projection_s01_vs_s02 = dataSigBESIII[i].makeProjection(Projection2D(data_projection_s01, data_projection_s02), PlotOptions::Prefix(("Data_s01_vs_s02_" + a[0]).c_str()));
         INFO("Done data Projections");
          auto mc_projection_s01 = QMI::proj_1D_psi3770(mcSig, mcTag, psi, s01_args_fit);
          auto mc_projection_s01_vs_s02 = QMI::proj_2D_psi3770(mcSig, mcTag, psi, s01_vs_s02_args_fit);
          TH1D* data_hist_s01 = dataSigBESIII[i].makeProjection(data_projection_s01, PlotOptions::Prefix(("Data_s01_" + a[0]).c_str()));
          TH1D* data_hist_s02 = dataSigBESIII[i].makeProjection(data_projection_s02, PlotOptions::Prefix(("Data_s02_" + a[0]).c_str()));

          data_hist_s01->Write();
          data_hist_s02->Write();
          data_projection_s01_vs_s02->Write();
          mc_projection_s01->Write();
          mc_projection_s01_vs_s02->Write();



        }
         int nFree = 0;
          for (auto& p : MPS){
              if (p->isFree()) nFree += 1;
          }

        TVectorD tvchi2(1);
          TVectorD tvndf(1);
          tvchi2[0] = _chi2;
          tvndf[0] = _nBins - nFree - 1;
          tvchi2.Write((a[0] + "_chi2").c_str());
          tvndf.Write((a[0] + "_ndf").c_str());


        chi2[i] = _chi2;
        nBins[i] = _nBins;
        INFO("_chi2 = "<<_chi2);
        INFO("_nBins = "<<_nBins);
        totalChi2 += _chi2;
        totalnBins += _nBins;
        i++;
      }


 
 
  fr.addChi2(totalChi2, totalnBins);  
    
  fr.print();
  fr.writeToFile(logFile);

  


  auto root_fr = new TFitResult(mini.fitResult());
  root_fr->SetName("FitResult");
  root_fr->Write();

  auto my_dd = [&A, &Abar](Event& evt){
    return QMI::dd(evt, A, Abar);
  };
  auto my_corr = [&compiledPoly, &MPS](Event& evt){
    return QMI::phaseCorrection(evt, compiledPoly, MPS);
  };

  QMI::writeValues(mcSig, my_dd, "dd");
  QMI::writeValues(mcSig, my_corr, "corr");

  QMI::writeDalitz(mcSig);



  tOutFile->Close();





  return 0;
}

