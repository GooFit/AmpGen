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
#include <fstream>
//#include <boost/algorithm/string/split.hpp>

using namespace AmpGen;

void readBinning(std::string binning, std::vector<real_t>& s01, std::vector<real_t>& s02, std::vector<int>& bins){
    std::ifstream inFile(binning.c_str());
    std::string line;
    while(std::getline(inFile, line)){
        std::vector<std::string> strVec;
        real_t x, y;
        int b;
        boost::algorithm::split(strVec, line, boost::is_any_of(" "));
        x = std::stod(strVec[0]);
        y = std::stod(strVec[1]);
        b = std::stoi(strVec[2]);
        s01.push_back(x);
        s02.push_back(y);
        bins.push_back(b);
    }
}

int nearestBinIdx(real_t s01, real_t s02, std::vector<real_t> s01_bin, std::vector<real_t> s02_bin){
  int idx = 0;
  real_t d2 = std::pow(s01_bin[idx] - s01, 2) + std::pow(s02_bin[idx] - s02, 2);
  
  for (int i=1;i<s01_bin.size();i++){

    real_t my_d2 = std::pow(s01_bin[i] - s01, 2) + std::pow(s02_bin[i] - s02, 2);
    if (my_d2 < d2) {
      idx = i;
      d2 = my_d2;
    }
  }
  return idx;
}


int nearestBinIdxOMP(real_t s01, real_t s02, std::vector<real_t> s01_bin, std::vector<real_t> s02_bin){
    struct Compare { float val; size_t index; };    
    #pragma omp declare reduction(minimum : struct Compare : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out)
    std::vector<real_t> d2(s01_bin.size());
    for (unsigned i=0;i<s01_bin.size();++i){
        d2[i] = std::pow(s01_bin[i] - s01, 2) +  std::pow(s02_bin[i] - s02, 2) ;
    }
    struct Compare min;
    min.val = d2[0];
    min.index = 0;
    #pragma omp parallel for reduction(minimum:min)
    for(int i = 1; i<d2.size(); i++) {
       if(d2[i]<min.val) { 
           min.val = d2[i];
           min.index = i;
        }
    }
    return min.index;
}


void binEvents(EventList list, std::vector<int>& binNum, std::vector<real_t> s01_bin, std::vector<real_t> s02_bin, std::vector<int> bins){
    for (auto& evt:list){
        real_t s01 = evt.s(0, 1);
        real_t s02 = evt.s(0, 2);
        int binIdx = nearestBinIdx(s01, s02, s01_bin, s02_bin);
        int bin = bins[binIdx];
        if (s02 < s01){
          bin = -bin;
        }
        binNum.push_back(bin);
    }
}

size_t countPerBin(std::vector<int> binnedEvents, int binNum){
    size_t r =0 ;
    for (int i=0;i<binnedEvents.size();++i){
        if (binnedEvents[i] == binNum){
            r++;
        }
    }
    return r;
}

size_t countPerBinDT(std::vector<int> binnedEvents, std::vector<int> binnedEvents_tag, int binNum, int binNum_tag ){
  size_t r =0;
  for (int i=0;i<binnedEvents.size();++i){
    if (binnedEvents[i]==binNum && binnedEvents_tag[i]==binNum_tag){
      r++;
    }
  }
  return r;
}

real_t expectCP(real_t F, real_t Fbar, real_t c, int CP){
  
    real_t n =  F + Fbar - 2 * CP * std::sqrt(F * Fbar)* c;
//    INFO("mu"<<CP<<"("<<c<<", "<<F<<", "<<Fbar<<") = "<<n);

    return n;
}
real_t expectDT(real_t F1, real_t F2, real_t Fbar1, real_t Fbar2, real_t c1, real_t c2, real_t s1, real_t s2){
    real_t n =  F1 * Fbar2 + F2 * Fbar1 - 2 * std::sqrt(F1 * F2 * Fbar1 * Fbar2) * (c1 * c2 + s1 * s2);
//    INFO("mu("<<c1<<", "<<", "<<s1<<", "<<F1<<", "<<Fbar1<<", "<<c2<<", "<<", "<<s2<<", "<<F2<<", "<<Fbar2<<") = "<<n);
    return n;
}
real_t expectCKM(real_t F, real_t Fbar, real_t c, real_t s, real_t x, real_t y, int sign){
  if (sign>0){
    return Fbar + (x * x + y * y) * F + 2 * std::sqrt(F * Fbar) * (c * x - s * y);
  }
  else{
    return F + (x * x + y * y) * Fbar + 2 * std::sqrt(F * Fbar) * (c * x + s * y);
  }



}
real_t chi2(int O, int E){
    return std::pow(((real_t)O - (real_t)E), 2)/(real_t)E;
}

real_t logPoisson(real_t x, real_t m){
  real_t y = 0;

//    y = -m + x * std::log(m) - std::lgamma(x + 1);
//  y = -0.5*std::pow(m - x, 2)/m - 0.5*std::log(M_PI * 2) - 0.5*std::log(m);
 real_t d = x;

 if (d==0) {
     y = 0;
     //y = - std::pow(m - x,2)/(1);///(m + x);
 }
 //y = -0.5* std::pow(m - x,2)/x;///(m + x);
 else{
 y = - std::pow(m - x,2)/(d);///(m + x);
 }
//    y = -m + x * std::log(m) - std::lgamma(x + 1);






 //y = - std::pow(m - x,2)/(x);///(m + x);
//  y = -0.5*std::pow(m - x,2)/(m);

 //y = -0.5* std::pow(m - x,2)/(x);///(m + x);
  //INFO("y = "<<y);
  return 0.5*y;
}
real_t logPoisson2(real_t x, real_t m, real_t err){
  real_t y = 0;

//    y = -m + x * std::log(m) - std::lgamma(x + 1);
//  y = -0.5*std::pow(m - x, 2)/m - 0.5*std::log(M_PI * 2) - 0.5*std::log(m);
 y = -0.5*std::pow( (m - x)/err,2);///(m + x);
//  y = -0.5*std::pow(m - x,2)/(m);

  //INFO("y = "<<y);
  return y;
}

real_t totalBinned(std::map<int, int> N){
  real_t norm =0;
  for (auto p : N){
    norm += p.second;
  }
  return norm;

}

std::map<int, real_t> normBinned(std::map<int, int> N){
  real_t norm =totalBinned(N);
  std::map<int, real_t> r;
  for (auto p:N){
    real_t f = (real_t)N[p.first]/norm;
    std::pair<int, real_t> pair({p.first, f});
    r.insert(pair);
  }
  return r;
}





int main(int argc, char * argv[]){
  OptionsParser::setArgs( argc, argv );
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  TRandom3 rndm(0);

  EventType sigType(NamedParameter<std::string>("EventType", "", "Signal Type to generate"));
  //size_t NInt(NamedParameter<size_t>("NInt", 1e7, "Number of events to calculate normalisation - should be large"));
  //size_t seed(NamedParameter<size_t>("Seed", 0, "Random seed for generation"));
//  size_t nEvents(NamedParameter<size_t>("nEvents", 1000, "number of events to generate, multiplies by the BR of each tag"));
  //size_t plot_nBins(NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms"));
  auto tags = NamedParameter<std::string>("TagTypes").getVector();
  auto btags = NamedParameter<std::string>("BTagTypes").getVector();
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  std::string logFile(NamedParameter<std::string>("LogFile", "Fit.log"));
  std::string psi3770Log(NamedParameter<std::string>("Psi3770Log", "Psi3770Fit.log"));
  std::string constantLog(NamedParameter<std::string>("constantLog", "constantFit.log"));
  std::string constrainedLog(NamedParameter<std::string>("constrainedLog", "constrainedFit.log"));
  std::string combinedLog(NamedParameter<std::string>("combinedLog", "combinedFit.log"));
//  std::string plotFile(NamedParameter<std::string>("Plots", "Fit.root"));
  #ifdef _OPENMP
    omp_set_num_threads( nThreads );
    INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
    omp_set_dynamic( 0 );
  #endif



  std::string BESIIIFile(NamedParameter<std::string>("BESIIIDataSample", "besiii.root", "Psi(3770) to DDbar events"));
  std::string LHCbFile(NamedParameter<std::string>("LHCbDataSample", "lhcb.root", "B+- to D h+- events"));

  std::string binning(NamedParameter<std::string>("binning", "KsPiPi_equal.txt"));
  std::vector<real_t> s01_bin, s02_bin;
  std::vector<int> bin;

  readBinning(binning, s01_bin, s02_bin, bin);
  std::map<std::string, std::map<int, int> > binnedEvents;
  std::map<std::pair<int, int>, int> binnedEventsDT;
  for(unsigned i=0;i<tags.size();++i ){
    auto a = split(tags[i], ' ');
    EventType tagType(Particle(a[1], {}, false).eventType());
    bool sameType = sigType == tagType;
  //  BESIIISameType[i] = sameType;
    EventList dataSig(BESIIIFile + ":Signal_" + a[0], sigType);
    EventList dataTag(BESIIIFile + ":Tag_" + a[0], tagType);
    std::vector<int> dataBins;
    std::map<int, int> binnedEvents_tag;
    binEvents(dataSig, dataBins, s01_bin, s02_bin, bin);

    if (sameType){
      std::vector<int> dataBins_tag;
      binEvents(dataTag, dataBins_tag, s01_bin, s02_bin, bin);
      for (int j=1;j<9;j++){
        for (int k=1;k<9;k++){
          std::pair<int, int> bin_pp({j, k});
          int N_pp = countPerBinDT(dataBins, dataBins_tag, j, k);
          std::pair<int, int> bin_pm({j, -k});
          int N_pm = countPerBinDT(dataBins, dataBins_tag, j, -k);
          std::pair<int, int> bin_mp({-j, k});
          int N_mp = countPerBinDT(dataBins, dataBins_tag, -j, k);
          std::pair<int, int> bin_mm({-j, -k});
          int N_mm = countPerBinDT(dataBins, dataBins_tag, -j, -k);
          std::pair<std::pair<int, int>, int> p_pp({bin_pp, N_pp});
          std::pair<std::pair<int, int>, int> p_pm({bin_pm, N_pm});
          std::pair<std::pair<int, int>, int> p_mp({bin_mp, N_mp});
          std::pair<std::pair<int, int>, int> p_mm({bin_mm, N_mm});
          binnedEventsDT.insert(p_pp);
          binnedEventsDT.insert(p_pm);
          binnedEventsDT.insert(p_mp);
          binnedEventsDT.insert(p_mm);
        }
      }
    }
    else{
      for (int j=1;j<9;j++){
        int Np = countPerBin(dataBins, j);
        int Nm = countPerBin(dataBins, -j);
        std::pair<int, int> p_p({j, Np});
        std::pair<int, int> p_m({-j, Nm});
        binnedEvents_tag.insert(p_p);
        binnedEvents_tag.insert(p_m);
      }
      binnedEvents.insert(std::pair<std::string, std::map<int, int> >({a[0], binnedEvents_tag}));
    }
  }
  for(unsigned i=0;i<btags.size();i++){
    auto a = split(btags[i], ' ');
    EventList data(LHCbFile + ":" + a[0], sigType);
    std::vector<int> dataBins;
    std::map<int, int> binnedEvents_tag;
    binEvents(data, dataBins, s01_bin, s02_bin, bin);
    for (int j=1;j<9;j++){
      int Np = countPerBin(dataBins, j);
      int Nm = countPerBin(dataBins, -j);
      std::pair<int, int> p_p({j, Np});
      std::pair<int, int> p_m({-j, Nm});
      binnedEvents_tag.insert(p_p);
      binnedEvents_tag.insert(p_m);
    }
    binnedEvents.insert(std::pair<std::string, std::map<int, int> >({a[0], binnedEvents_tag}));
  }

  for (auto p : binnedEvents){
    for (auto q : p.second){
      INFO(p.first<<" "<<q.first<< " = "<<q.second);
    }
  }

  auto F = normBinned(binnedEvents["Kppim"]);
  auto Fbar = normBinned(binnedEvents["Kmpip"]);
  real_t totalKK = totalBinned(binnedEvents["KK"]);
  real_t totalKppim = totalBinned(binnedEvents["Kppim"]);
  real_t totalKmpip = totalBinned(binnedEvents["Kmpip"]);
  real_t totalKspi0 = totalBinned(binnedEvents["Kspi0"]);

  for (int i=1;i<9;++i){
      //MPS["F"+std::to_string(i)]->setCurrentFitVal(F[i]);
      MPS["F"+std::to_string(i)]->setCurrentFitVal( (binnedEvents["Kppim"][i] + binnedEvents["Kmpip"][-i])/(totalKppim + totalKmpip));   ;
      //MPS["F"+std::to_string(-i)]->setCurrentFitVal(Fbar[i]);
      //MPS["F"+std::to_string(-i)]->setCurrentFitVal(Fbar[i]);
      MPS["F"+std::to_string(-i)]->setCurrentFitVal( (binnedEvents["Kppim"][-i] + binnedEvents["Kmpip"][i])/(totalKppim + totalKmpip));   ;
  }

//  real_t totalKspipi = totalBinned(binnedEventsDT);

for (auto p : F){
  INFO("F"<<p.first<<" = "<<p.second);
}

for (auto p : Fbar){
  INFO("Fbar"<<p.first<<" = "<<p.second);
}
auto min_Kppim = [&binnedEvents, &MPS, &totalKppim, &totalKmpip](){
        std::map<int, real_t> mu;
        real_t norm =0;
        for (int i=1;i<9;++i){
            real_t F_i = MPS["F" + std::to_string(i)]->mean();
            real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();
            real_t mu_i = F_i;
            real_t mu_mi = Fbar_i;
            mu.insert(std::pair<int, real_t>({i, mu_i}));
            mu.insert(std::pair<int, real_t>({-i, mu_mi}));
            norm += mu_i + mu_mi;
        }
        real_t ll =0;
        for (int i=1;i<9;++i){
          real_t F_i = MPS["F" + std::to_string(i)]->mean();
            real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();

            //real_t m = (mu[i] + mu[-1])/norm;
///            real_t m = F_i;
            real_t E = totalKppim * F_i + totalKppim * Fbar_i;
            real_t n = binnedEvents["Kppim"][i] + binnedEvents["Kppim"][-i];
            ll += logPoisson(n, E);
//            m = Fbar_i;
//            E = totalKppim * m;
//            n = binnedEvents["Kppim"][-i];
//            ll += logPoisson(n, E);

        }
        return -2 * ll;
    };
    auto min_Kmpip = [&binnedEvents, &MPS, &totalKmpip, &totalKppim](){
        std::map<int, real_t> mu;
        real_t norm =0;
        for (int i=1;i<9;++i){
            real_t F_i = MPS["F" + std::to_string(i)]->mean();
            real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();
            real_t mu_i = Fbar_i;
            real_t mu_mi = F_i;
            mu.insert(std::pair<int, real_t>({i, mu_i}));
            mu.insert(std::pair<int, real_t>({-i, mu_mi}));
            norm += mu_i + mu_mi;
        }
        real_t ll =0;
        for (int i=1;i<9;++i){
              real_t F_i = MPS["F" + std::to_string(i)]->mean();
            real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();


            //real_t m = (mu[i] + mu[-1])/norm;
 //           real_t m = Fbar_i;
            real_t E = totalKmpip * Fbar_i + totalKmpip * F_i;
            //real_t n = binnedEvents["Kmpip"][i] + binnedEvents["Kmpip"][-i];
            real_t n = binnedEvents["Kmpip"][i] + binnedEvents["Kmpip"][-i];
            ll += logPoisson(n, E);
//            m = F_i;
//            E = totalKmpip* m;
//            n = binnedEvents["Kmpip"][-i];
//            ll += logPoisson(n, E);
        }
        return -2 * ll;
    };


  auto min_KK = [&binnedEvents, &MPS, &F, &Fbar, &totalKK](){
    std::map<int, real_t> mu;
    real_t norm = 0;
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
    real_t F_i = MPS["F" + std::to_string(i)]->mean();
    real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();


      //real_t mu_i = expectCP(F[i], Fbar[i], ci, 1);
      real_t mu_i = expectCP(F_i, Fbar_i, ci, 1);


      norm += mu_i;
      mu.insert(std::pair<int, real_t>({i, mu_i}));
    //real_t mu_mi = expectCP(F[-i], Fbar[-i], ci, 1);
    real_t mu_mi = expectCP(Fbar_i, F_i, ci, 1);
    norm += mu_mi;
     // int exp_mi = mu_mi * totalKK;
      mu.insert(std::pair<int, real_t>({-i, mu_mi}));
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t dci = MPS["c" + std::to_string(i)]->err();
      real_t m = (mu[i] + mu[-i])/norm;
      real_t E = totalKK * m; 
      real_t n = binnedEvents["KK"][i] + binnedEvents["KK"][-i];
      //ll += logPoisson2(n, E, std::sqrt(n + 1));
      ll += logPoisson(n, E);
      
      real_t exp_p = totalKK * (F[i] + Fbar[i] - 2 * ci * std::sqrt(F[i] * Fbar[i]));
      real_t exp_m = totalKK * (F[-i] + Fbar[-i] - 2 * ci * std::sqrt(F[-i] * Fbar[-i]));
      
    }
    return -2 * ll;
  };
  auto min_Kspi0 = [&binnedEvents, &MPS, &F, &Fbar, &totalKspi0](){
    std::map<int, real_t> mu;
    real_t norm = 0;
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
    real_t F_i = MPS["F" + std::to_string(i)]->mean();
    real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();


      //real_t mu_i = expectCP(F[i], Fbar[i], ci, -1);
      real_t mu_i = expectCP(F_i, Fbar_i, ci, -1);

      mu.insert(std::pair<int, real_t>({i, mu_i}));
      norm += mu_i;
    //real_t mu_mi = expectCP(F[-i], Fbar[-i], ci, -1);
    real_t mu_mi = expectCP(Fbar_i, F_i, ci, -1);

      mu.insert(std::pair<int, real_t>({-i, mu_mi}));
      norm += mu_mi;
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){

      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t dci = MPS["c" + std::to_string(i)]->err();
      real_t m = (mu[i] + mu[-i])/norm;
      real_t E = m * totalKspi0;
      real_t n = binnedEvents["Kspi0"][i] + binnedEvents["Kspi0"][-i];
      
      //ll += logPoisson2(n, E, std::sqrt(n + 1));
      ll += logPoisson(n, E);
      
    }
    return -2 * ll;
  };


  real_t totalKspipi = 0;
  for (auto p : binnedEventsDT){
    totalKspipi += p.second;
  }

  auto min_Kspipi = [&binnedEventsDT, &MPS, &F, &Fbar, &totalKspipi](){
    std::map<std::pair<int, int>, real_t> mu;
    real_t norm =0;
    for (int i=1;i<9;i++){
    real_t F_i = MPS["F" + std::to_string(i)]->mean();
    real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();


      real_t ci = MPS["c" + std::to_string(i)]->mean();

      real_t si = MPS["s" + std::to_string(i)]->mean();

      for (int j=1;j<9;j++){
    real_t F_j = MPS["F" + std::to_string(j)]->mean();
    real_t Fbar_j = MPS["F" + std::to_string(-j)]->mean();


        real_t cj = MPS["c" + std::to_string(j)]->mean();

        real_t sj = MPS["s" + std::to_string(j)]->mean();


        std::pair<int, int> b_pp({i, j});
        std::pair<int, int> b_pm({i, -j});
        std::pair<int, int> b_mp({-i, j});
        std::pair<int, int> b_mm({-i, -j});
//        real_t mu_pp = expectDT(F[i], F[j], Fbar[i], Fbar[j], ci, cj, si, sj);
//        real_t mu_pm = expectDT(F[i], F[-j], Fbar[i], Fbar[-j], ci, cj, si, -sj);
//        real_t mu_mp = expectDT(F[-i], F[j], Fbar[-i], Fbar[j], ci, cj, -si, sj);
//        real_t mu_mm = expectDT(F[-i], F[-j], Fbar[-i], Fbar[-j], ci, cj, -si, -sj);
        real_t mu_pp = expectDT(F_i, F_j, Fbar_i, Fbar_j, ci, cj, si, sj);
        real_t mu_pm = expectDT(F_i, Fbar_j, Fbar_i, F_j, ci, cj, si, -sj);
        real_t mu_mp = expectDT(Fbar_i, F_j, F_i, Fbar_j, ci, cj, -si, sj);
        real_t mu_mm = expectDT(Fbar_i, Fbar_j, F_i, F_j, ci, cj, -si, -sj);




        std::pair<std::pair<int, int>, real_t> p_pp({b_pp, mu_pp});
        std::pair<std::pair<int, int>, real_t> p_pm({b_pm, mu_pm});
        std::pair<std::pair<int, int>, real_t> p_mp({b_mp, mu_mp});
        std::pair<std::pair<int, int>, real_t> p_mm({b_mm, mu_mm});

        mu.insert(p_pp);
        mu.insert(p_pm);
        mu.insert(p_mp);
        mu.insert(p_mm);
        norm +=  mu_pp + mu_pm + mu_mp + mu_mm;
      }
    }
    real_t ll = 0;
    for (int i=1;i<9;i++){
      for (int j=1;j<9;j++){
        std::pair<int, int> b_pp({i, j});
        std::pair<int, int> b_pm({i, -j});
        std::pair<int, int> b_mp({-i, j});
        std::pair<int, int> b_mm({-i, -j});
        real_t n_p = binnedEventsDT[b_pp] + binnedEventsDT[b_mm];
        real_t n_m = binnedEventsDT[b_pm] + binnedEventsDT[b_mp];
        real_t mu_p = (mu[b_pp] + mu[b_mm])/norm;
        real_t mu_m = (mu[b_pm] + mu[b_mp])/norm;
        real_t E_p = mu_p * totalKspipi;
        real_t E_m = mu_m * totalKspipi;
        real_t dci = MPS["c" + std::to_string(i)]->err();
        real_t dcj = MPS["c" + std::to_string(j)]->err();
        real_t dsi = MPS["s" + std::to_string(i)]->err();
        real_t dsj = MPS["s" + std::to_string(j)]->err();
        real_t err = std::sqrt(dci * dci + dcj * dcj + dsi * dsi + dsj * dsj);
        //ll += logPoisson2(n_p, E_p, std::sqrt(E_p + 1));
//        ll += logPoisson(n_p + n_m, E_p + E_m);
        ll += logPoisson(n_p, E_p);
        ll += logPoisson(n_m, E_m);
        //ll += logPoisson2(n_m, E_m, std::sqrt(E_m + 1));
        
      }
    }
    return -2 * ll;
  };
  auto min_CP = [&min_KK, &min_Kspi0](){//, &min_Kppim, &min_Kmpip](){
    return min_KK() + min_Kspi0(); //+ min_Kppim() + min_Kmpip();
  };
  auto min_psi3770 = [&min_CP, &min_Kspipi](){
    return min_CP() + min_Kspipi();
  };

  MinuitParameterSet * MPS_psi3770 = new MinuitParameterSet();
  for (int i=1;i<9;i++){
    MPS_psi3770->add(MPS["c" + std::to_string(i)]);
    MPS_psi3770->add(MPS["s" + std::to_string(i)]);
    MPS_psi3770->add(MPS["F" + std::to_string(i)]);
    MPS_psi3770->add(MPS["F" + std::to_string(-i)]);
  }

  Minimiser mini_KK(min_KK, MPS_psi3770);
  Minimiser mini_Kspi0(min_Kspi0, MPS_psi3770);
  Minimiser mini_CP(min_CP, MPS_psi3770);

  Minimiser mini_psi3770(min_psi3770, MPS_psi3770);
  real_t ll_KK = min_KK();
  real_t ll_Kspi0 = min_Kspi0();
//  mini_CP.gradientTest();
//  mini_KK.gradientTest();
//  mini_KK.doFit();
//  mini_Kspi0.doFit();
// mini_CP.doFit();

//  MPS["CKM::x+"]->fix();
//  MPS["CKM::y+"]->fix();
//  MPS["CKM::x-"]->fix();
//  MPS["CKM::y-"]->fix();
  mini_psi3770.gradientTest();
  mini_psi3770.doFit();
  FitResult * fr = new FitResult(mini_psi3770);
  fr->writeToFile(psi3770Log);
//  delete MPS_psi3770; 
 // return 0;
//  if (psi3770Log != ""){
   INFO("Getting mu/sigma for parameters");
   auto meanAndErr_psi3770 = QMI::fitValAndErr(psi3770Log);
   INFO("Getting inverse Covariance matrix");
   auto invCovMatrix_psi3770 = QMI::invCovarianceMatrix(psi3770Log);
   INFO("Done invConv");
   auto constraint = [&meanAndErr_psi3770, &invCovMatrix_psi3770, &MPS](){
     return QMI::myGaussConstraint(meanAndErr_psi3770, invCovMatrix_psi3770, MPS);
   };
   INFO("Built Constraint as function");
   //}
  
  real_t totalBp2DKp = totalBinned(binnedEvents["Bp2DKp"]);
  auto min_Bp2DKp = [&binnedEvents, &MPS, &F, &Fbar, &totalBp2DKp](){
    std::map<int, real_t> mu;
    real_t norm = 0;
    real_t x = MPS["CKM::x+"]->mean();
    real_t y = MPS["CKM::y+"]->mean();
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t si = MPS["s" + std::to_string(i)]->mean();
      //real_t mu_i = expectCP(F[i], Fbar[i], ci, -1);
      real_t mu_i = expectCKM(F[i], Fbar[i], ci, si, x, y, 1);

      mu.insert(std::pair<int, real_t>({i, mu_i}));
      norm += mu_i;
    real_t mu_mi = expectCKM(F[-i], Fbar[-i], ci, -si, x, y, 1);

      mu.insert(std::pair<int, real_t>({-i, mu_mi}));
      norm += mu_mi;
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){

      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t dci = MPS["c" + std::to_string(i)]->err();
      real_t m = (mu[i])/norm;
      real_t E = m * totalBp2DKp;
      real_t n = binnedEvents["Bp2DKp"][i];// + binnedEvents["Bp2DKp"][-i];
      
      //ll += logPoisson2(n, E, std::sqrt(n + 1));
      

      real_t mm = mu[-i]/norm;
      real_t Em = mm * totalBp2DKp;
      real_t nm = binnedEvents["Bp2DKp"][-i];
      ll += logPoisson(n, E);
      ll += logPoisson(nm, Em);
 //     ll += logPoisson(n+ nm, E + Em);
      
    }
    return -2 * ll;
  };
  real_t totalBm2DKm = totalBinned(binnedEvents["Bm2DKm"]);
  auto min_Bm2DKm = [&binnedEvents, &MPS, &F, &Fbar, &totalBm2DKm](){
    std::map<int, real_t> mu;
    real_t norm = 0;
    real_t x = MPS["CKM::x-"]->mean();
    real_t y = MPS["CKM::y-"]->mean();
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t si = MPS["s" + std::to_string(i)]->mean();
      //real_t mu_i = expectCP(F[i], Fbar[i], ci, -1);
      real_t mu_i = expectCKM(F[i], Fbar[i], ci, si, x, y, -1);

      mu.insert(std::pair<int, real_t>({i, mu_i}));
      norm += mu_i;
    real_t mu_mi = expectCKM(F[-i], Fbar[-i], ci, -si, x, y, -1);

      mu.insert(std::pair<int, real_t>({-i, mu_mi}));
      norm += mu_mi;
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){

      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t dci = MPS["c" + std::to_string(i)]->err();
      real_t m = (mu[i])/norm;
      real_t E = m * totalBm2DKm;
      real_t n = binnedEvents["Bm2DKm"][i]; //+ binnedEvents["Bm2DKm"][-i];
      
      //ll += logPoisson2(n, E, std::sqrt(n + 1));

      real_t mm = mu[-i]/norm;
      real_t Em = mm * totalBm2DKm;
      real_t nm = binnedEvents["Bm2DKm"][-i];
      ll += logPoisson(n, E);
      ll += logPoisson(nm, Em);
//      ll += logPoisson(n + nm, E + Em);
 
      
    }
    return -2 * ll;
  };




  for (int i=1;i<9;i++){
    MPS["c" + std::to_string(i)]->fix();
    MPS["s" + std::to_string(i)]->fix();
  }
  MPS["CKM::x+"]->setFree();
  MPS["CKM::y+"]->setFree();
  MPS["CKM::x-"]->setFree();
  MPS["CKM::y-"]->setFree();


  auto min_CKM = [&min_Bp2DKp, &min_Bm2DKm](){
    return min_Bp2DKp() + min_Bm2DKm();
  };



  Minimiser mini_CKM(min_CKM, &MPS);
  mini_CKM.doFit();
  FitResult * fr_CKM = new FitResult(mini_CKM);
  fr_CKM->writeToFile(constantLog);
 


  for (int i=1;i<9;i++){
    MPS["c" + std::to_string(i)]->setFree();
    MPS["s" + std::to_string(i)]->setFree();
  }

  auto min_CKM_constrained = [&min_CKM, &constraint](){
    return min_CKM() + constraint();
  };
  Minimiser mini_constrained(min_CKM_constrained, &MPS);
  mini_constrained.doFit();
  FitResult * fr_constrained = new FitResult(mini_constrained);
  fr_constrained->writeToFile(constrainedLog);
 

  auto min_comb = [&min_psi3770, &min_CKM](){
    return min_psi3770() + min_CKM();
  };
 

  Minimiser mini_comb(min_comb, &MPS);
  mini_comb.doFit();
  FitResult * fr_combined = new FitResult(mini_comb);
  fr_combined->writeToFile(combinedLog);
  
    return 0;
}
