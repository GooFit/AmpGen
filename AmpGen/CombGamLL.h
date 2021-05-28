#ifndef COMBGAMLL
#define COMBGAMLL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/Generator.h"


#include <tuple>
#include <TRandom3.h>
namespace AmpGen
{
  class EventList;
  /**
   * @class CombCorrLL
   * @brief A combined log-likelihood for correlated amplitudes.
   **/
  class CombGamLL
  {
      private:
        std::vector<EventList> m_SigData;
        EventType m_SigType;
        EventList m_MC;
        MinuitParameterSet m_MPS;
        std::vector<pCoherentSum> m_Psi;
    
      public:

        CombGamLL(EventType SigType, std::vector<EventList> SigData, MinuitParameterSet MPS, std::vector<int> gammaSign,
         std::vector<bool> useXY, std::vector<bool> Conj, size_t seed, size_t NInt ):
            m_SigData(SigData),
            m_SigType(SigType),
            m_MPS(MPS) {
                TRandom3 rndm;
                rndm.SetSeed( seed );
                gRandom = &rndm;

                
                m_MC =  Generator<>(SigType, &rndm).generate(NInt);
                for (size_t i=0; i < m_SigData.size(); i++){

                    auto psi = pCoherentSum(SigType, MPS, gammaSign[i], useXY[i], Conj[i]);
                    psi.setEvents(m_SigData[i]);
                    psi.setMC(m_MC);
                    psi.prepare();
                    m_Psi.push_back(psi);
                }
            }

        CombGamLL(EventType SigType, std::vector<EventList> SigData, MinuitParameterSet MPS, std::vector<int> gammaSign,
         std::vector<bool> useXY, std::vector<bool> Conj, EventList SigInt ):
            m_SigData(SigData),
            m_SigType(SigType),
            m_MC(SigInt),
            m_MPS(MPS) {
            for (size_t i=0; i < m_SigData.size(); i++){

                    auto psi = pCoherentSum(SigType, MPS, gammaSign[i], useXY[i], Conj[i]);
                    psi.setEvents(m_SigData[i]);
                    psi.setMC(m_MC);
                    psi.prepare();
                    m_Psi.push_back(psi);
                }
            }

        void reset(){
            for (size_t i=0;i<m_Psi.size();i++) m_Psi[i].reset(true);
        }

        real_t LL(int i){
            real_t n = m_Psi[i].norm();
            real_t ll =0;
            for (size_t j=0;j<m_SigData[i].size();j++){
                real_t pdf = std::norm(m_Psi[i].getVal(m_SigData[i][j]))/n;
                ll += -2 * log(pdf);
            }
            return ll;
        }

        real_t getVal(){
            real_t ll=0;
            for (size_t i=0;i<m_SigData.size();i++){
                ll+=LL(i);
            }
            return ll;
        }
        void makeProjection(size_t i, std::string plotFile, std::string prefix, size_t nBins){
            TFile * f = TFile::Open(plotFile.c_str(), "UPDATE");
            real_t norm = m_Psi[i].norm();
            auto projections = m_SigType.defaultProjections(nBins);
            for (auto& proj : projections){
                auto data = proj(m_SigData[i], Prefix(prefix));
                auto hist = proj.plot(prefix);
                real_t integral = 0;
                for (size_t j=0;j<m_MC.size();j++){
                    real_t p =  std::norm(m_Psi[i].getVal(m_MC[j]))/norm ;
                    hist->Fill( proj(m_MC[j]), p);
                    integral += p;
                }
                hist->Scale(data->Integral()/integral);
                hist->SetName( (std::string("MC_") + hist->GetName()).c_str());
                hist->Write();
                data->Write();


                auto bins = hist->GetBin(hist->GetEntries()) -1;
                auto x0 = hist->GetBinCenter(0);
                auto x1 = hist->GetBinCenter(bins +1);

                TH1D * pull = new TH1D( (std::string("Pull_") + hist->GetName()).c_str(), "Pull", bins, x0, x1 );
                for (int j=0;j<hist->GetEntries();j++){
                    double p=0;
                    double d = hist->GetBinContent(j) - data->GetBinContent(j);
                    double s2 = hist->GetBinContent(j) + data->GetBinContent(j);
                    double x = hist->GetBinCenter(j);
                    if (s2!=0) p=d/sqrt(s2);
                    pull->Fill(x, p);
                }
                pull->Write();
            }
            for (size_t j=0; j < projections.size(); j++){
                for(size_t k=j+1; k < projections.size(); k++){
                    auto dalitz = Projection2D(projections[j], projections[k]);
                    auto hdalitz = dalitz.plot(prefix);
                    auto data = m_SigData[i].makeProjection(dalitz, Prefix(prefix));
                    real_t integral=0;
                    for (size_t l=0;l<m_MC.size();l++){

                        real_t p =  std::norm(m_Psi[i].getVal(m_MC[l]))/norm ;
                        auto pos = dalitz(m_MC[l]);
                        hdalitz->Fill(pos.first, pos.second,p);
                        integral += p;
                    }
                    hdalitz->Scale(data->Integral() / integral);
                    hdalitz->SetName( (std::string("MC_") + hdalitz->GetName()).c_str() );
                    hdalitz->Write();
                    data->Write();

                    auto px = hdalitz->ProjectionX();
                    auto py = hdalitz->ProjectionY();

                    auto bins = px->GetBin(px->GetEntries() -1);
                    auto x0 = px->GetBinCenter(0);
                    auto x1 = px->GetBinCenter(bins +1);

                    bins = py->GetBin(py->GetEntries() -1);
                    auto y0 = py->GetBinCenter(0);
                    auto y1 = py->GetBinCenter(bins +1);


                 //   bins = int(sqrt(bins));


                    TH2D * pull_2D = new TH2D( (std::string("Pull_") + hdalitz->GetName() ).c_str(), "Pull", bins, x0, x1, bins, y0, y1 );
                    for (int l=0;l<data->GetEntries();l++){
                        for (int m=0;m<hdalitz->GetEntries();m++){
                            double p=0;
                            double d = hdalitz->GetBinContent(l,m) - data->GetBinContent(l,m);
                            double s2 = hdalitz->GetBinContent(l,m) + data->GetBinContent(l,m);
                            if (s2!=0) p = d/sqrt(s2);
                            double x = px->GetBinCenter(l);
                            double y = py->GetBinCenter(m);
                            pull_2D->Fill(x, y, p);
                        }
                    }
                    pull_2D->Write();

                    
                }
            }
            f->Close();
        }
            
        };
//              auto psi = pCoherentSum(SigType[0], MPS, BgammaSigns[i], BuseXYs[i], BConj[i]); 

}


#endif