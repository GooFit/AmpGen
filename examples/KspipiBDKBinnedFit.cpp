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
    return F + Fbar - 2 * CP * std::pow(F * Fbar, 0.5)* c;
}
real_t expectDT(real_t F1, real_t F2, real_t Fbar1, real_t Fbar2, real_t c1, real_t c2, real_t s1, real_t s2){
    return F1 * Fbar2 + F2 * Fbar1 - 2 * std::pow(F1 * F2 * Fbar1 * Fbar2, 0.5) * (c1 * c2 - s1 * s2);
}

real_t expectCKM(real_t F, real_t Fbar, real_t c, real_t s, real_t x, real_t y, int sign){
  if (sign>0){
    return Fbar + (x * x + y * y) * F + 2 * std::sqrt(F * Fbar) * (c * x - s * y);
  }
  else{
    return F + (x * x + y * y) * Fbar + 2 * std::sqrt(F * Fbar) * (c * x + s * y);
  }



}

real_t chi2(real_t O, real_t E){
    return std::pow((O - E), 2)/E;
}

real_t logPoisson2(size_t x, size_t m){
  return -m + m * std::log(x) - std::lgamma(x + 1);
}


real_t logPoisson(real_t x, real_t m){
  real_t y = 0;

 
//    y = -m + x * std::log(m) - std::lgamma(x + 1);
//  y = -0.5*std::pow(m - x, 2)/m - 0.5*std::log(M_PI * 2) - 0.5*std::log(m);
 if (x != 0){
 y = -std::pow(m - x,2)/x;///(m + x);
 
//    y = -m + x * std::log(m) - std::lgamma(x + 1);
 }
 else{

 y = -std::pow(m - x,2);///(m + x);
 }
//  y = -0.5*std::pow(m - x,2)/(m);
 //y = -std::pow(m - x,2)/x;///(m + x);
  //INFO("y = "<<y);
  return 0.5*y;
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
  //auto tags = NamedParameter<std::string>("TagTypes").getVector();
  auto btags = NamedParameter<std::string>("BTagTypes").getVector();
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  std::string logFile(NamedParameter<std::string>("LogFile", "Fit.log"));
  std::string plotFile(NamedParameter<std::string>("Plots", "Fit.root"));
  #ifdef _OPENMP
    omp_set_num_threads( nThreads );
    INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
    omp_set_dynamic( 0 );
  #endif


  std::string psi3770Log(NamedParameter<std::string>("psi3770Log", ""));

//  std::string BESIIIFile(NamedParameter<std::string>("BESIIIDataSample", "besiii.root", "Psi(3770) to DDbar events"));
  std::string LHCbFile(NamedParameter<std::string>("LHCbDataSample", "lhcb.root", "B+- to D h+- events"));

  std::string binning(NamedParameter<std::string>("binning", "KsPiPi_equal.txt"));
  std::vector<real_t> s01_bin, s02_bin;
  std::vector<int> bin;

  readBinning(binning, s01_bin, s02_bin, bin);
  std::map<std::string, std::map<int, int> > binnedEvents;
  //std::map<std::pair<int, int>, int> binnedEventsDT;
  
//  auto F = normBinned(binnedEvents["Kppim"]);
//  auto Fbar = normBinned(binnedEvents["Kmpip"]);
//  real_t totalKK = totalBinned(binnedEvents["KK"]);
//  real_t totalKspi0 = totalBinned(binnedEvents["Kspi0"]);
//  real_t totalKspipi = totalBinned(binnedEventsDT);
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
  std::map<int,real_t> F;
  std::map<int,real_t> Fbar;
  for (int i=1;i<9;i++){
    
    std::pair<int, real_t> p_p({i, MPS["F" + std::to_string(i)]->mean()});
    std::pair<int, real_t> p_m({-i, MPS["F" + std::to_string(-i)]->mean()});
    std::pair<int, real_t> pbar_p({i, MPS["Fbar" + std::to_string(i)]->mean()});
    std::pair<int, real_t> pbar_m({-i, MPS["Fbar" + std::to_string(-i)]->mean()});
    F.insert(p_p);
    F.insert(p_m);
    Fbar.insert(pbar_p);
    Fbar.insert(pbar_m);
  }

  real_t totalBp = 0;
  real_t totalBm = 0;
  for (auto p : binnedEvents["Bp2DKp"]){
    totalBp += p.second;
  }
  for (auto p : binnedEvents["Bm2DKm"]){
    totalBm += p.second;
  }

  for (auto p : binnedEvents){
    for (auto q : p.second){
      INFO(p.first<<" "<<q.first<< " = "<<q.second);
    }
  }

   auto min_Bp = [&binnedEvents, &MPS, &F, &Fbar, &totalBp](){
    std::map<int, real_t> mu;
    real_t x = MPS["CKM::x+"]->mean();
    real_t y = MPS["CKM::y+"]->mean();
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t si = MPS["s" + std::to_string(i)]->mean();
      real_t Fi = F[i];
      real_t Fmi = F[-i];
      real_t Fbari = Fbar[i];
      real_t Fbarmi = Fbar[-i];
      real_t F_i = MPS["F" + std::to_string(i)]->mean();
      real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();



      //real_t mu_i = expectCKM(Fi, Fbari, ci, si, x, y, 1);
      real_t mu_i = expectCKM(F_i, Fbar_i, ci, si, x, y, 1);
      real_t exp_i = mu_i * totalBp;
      mu.insert(std::pair<int, real_t>({i, exp_i}));
      //real_t mu_mi = expectCKM(Fmi, Fbarmi, ci, -si, x, y, 1);
      real_t mu_mi = expectCKM(Fbar_i, F_i, ci, -si, x, y, 1);
      real_t exp_mi = mu_mi * totalBp;
      mu.insert(std::pair<int, real_t>({-i, exp_mi}));
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){
   //   real_t m = mu[i] + mu[-i];
//      int n = binnedEvents["Kspi0"][i] + binnedEvents["Kspi0"][-i];
      ll += logPoisson(binnedEvents["Bp2DKp"][i], mu[i]);
      ll += logPoisson(binnedEvents["Bp2DKp"][-i], mu[-i]);
      //ll += logPoisson(binnedEvents["Bp2DKp"][i] + binnedEvents["Bp2DKp"][-i], mu[i] + mu[-i]);
    }
    return  ll;
  };
  auto min_Bm = [&binnedEvents, &MPS, &F, &Fbar, &totalBm](){
    std::map<int, real_t> mu;
    real_t x = MPS["CKM::x-"]->mean();
    real_t y = MPS["CKM::y-"]->mean();
    for (int i=1;i<9;++i){
      real_t ci = MPS["c" + std::to_string(i)]->mean();
      real_t si = MPS["s" + std::to_string(i)]->mean();
      real_t Fi = F[i]; 
      real_t Fmi = F[-i];
      real_t Fbari = Fbar[i];
      real_t Fbarmi = Fbar[-i]; 
      real_t F_i = MPS["F" + std::to_string(i)]->mean();
      real_t Fbar_i = MPS["F" + std::to_string(-i)]->mean();



      //real_t mu_i = expectCKM(Fi, Fbari, ci, si, x, y, -1);
      real_t mu_i = expectCKM(F_i, Fbar_i, ci, si, x, y, -1);
      real_t exp_i = mu_i * totalBm;
      mu.insert(std::pair<int, real_t>({i, exp_i}));
      //real_t mu_mi = expectCKM(Fmi, Fbarmi, ci, -si, x, y, -1);
      real_t mu_mi = expectCKM(Fbar_i, F_i, ci, -si, x, y, -1);
      real_t exp_mi = mu_mi * totalBm;
      mu.insert(std::pair<int, real_t>({-i, exp_mi}));
    }
    real_t ll =0 ;
    for (int i=1;i<9;++i){
     // int m = mu[i] + mu[-i];
//      int n = binnedEvents["Kspi0"][i] + binnedEvents["Kspi0"][-i];
      ll += logPoisson(binnedEvents["Bm2DKm"][i], mu[i]);
      ll += logPoisson(binnedEvents["Bm2DKm"][-i], mu[-i]);
//      ll += logPoisson(binnedEvents["Bm2DKm"][i] + binnedEvents["Bm2DKm"][-i], mu[i] + mu[-i]);
    }
    return  ll;
  };

  auto min_B = [&min_Bp, &min_Bm](){
    return min_Bp() + min_Bm();
  };




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
  auto min_B_constrained = [&min_B, &constraint](){
    return min_B() + constraint();
  };
  Minimiser mini (min_B_constrained, &MPS);
  mini.doFit();
 }
 else{
  Minimiser mini (min_B, &MPS);
  mini.doFit();
 }


 
  


  
    return 0;
}
