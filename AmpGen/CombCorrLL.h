#ifndef COMBCORRLL
#define COMBCORRLL

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/Generator.h"

#include <tuple>
namespace AmpGen
{
  class EventList;
  /**
   * @class CombCorrLL
   * @brief A combined log-likelihood for correlated amplitudes.
   **/
  class CombCorrLL
  {
      private:
        std::vector<EventList> m_SigData;
        std::vector<EventList> m_TagData;
        std::vector<EventList> m_SigInt;
        std::vector<EventList> m_TagInt;
        std::vector<EventType> m_SigType;
        std::vector<EventType> m_TagType;
        MinuitParameterSet m_mps;
        bool m_debug;
        std::vector<pCorrelatedSum> m_Psi = {};
        std::vector<pCorrelatedSum*> m_Psi2;
        std::vector<std::string> m_SumFactors;
        
        EventList m_MC;

      public:
        CombCorrLL() = default;
        CombCorrLL(std::vector<pCorrelatedSum>psi):
          m_Psi(psi)
          {
            for(int i=0;i<psi.size();i++){
              INFO("psi["<<i<<"]  ");//<<psi[i]);

             
              
            }

          }

//void makeProjection(pCorrelatedSum cs_tag, std::string tag_plotName, std::string prefix ,int nBins, EventList sigevents_tag, EventList tagevents_Tag, EventList sigMCevents_tag, EventList tagMCevents_tag){
void makeProjection(int i, std::string tag_plotName, std::string prefix ,int nBins){
        INFO("@"<<i);
        auto cs_tag = m_Psi[i];
        auto sigevents_tag = m_SigData[i];
        //auto sigMCevents_tag = m_SigInt[i];
        auto sigMCevents_tag = m_SigInt[0];
        auto tagevents_tag = m_TagData[i];
        auto tagMCevents_tag = m_TagInt[i];

        INFO( "norm[1] = " << cs_tag.norm() );
        TFile* f = TFile::Open(tag_plotName.c_str(),"UPDATE");
        auto projections = sigevents_tag.eventType().defaultProjections(nBins);
        real_t norm = cs_tag.norm();
        for( auto& projection : projections ){
          auto data_plot = projection(sigevents_tag, Prefix(prefix));
          auto hist = projection.plot(prefix);
          real_t integral = 0;
          for(unsigned i = 0 ; i != sigMCevents_tag.size(); ++i)
          {
            hist->Fill( projection( sigMCevents_tag[i] ), std::norm(cs_tag.getVal( sigMCevents_tag[i], tagMCevents_tag[i] ))/norm );
            integral += std::norm(cs_tag.getVal(sigMCevents_tag[i], tagMCevents_tag[i]))/norm;

          }
          auto data_integral = data_plot->Integral();
          auto hist_integral = hist->Integral();
          INFO("data integral = "<<data_integral);
          INFO("hist integral = "<<hist_integral);
          INFO("hist integral = "<<integral);
          //hist->Scale( data_plot->Integral() / hist->Integral() );
          hist->Scale( data_plot->Integral() / integral );
          hist->SetName( (std::string("MC_")+hist->GetName()).c_str() );
          hist->Write();
          data_plot->Write();

          auto bins = hist->GetBin(hist->GetEntries()) -1;
          auto x0 = hist->GetBinCenter(0);
          auto x1 = hist->GetBinCenter(bins +1);

          TH1D * pull = new TH1D( (std::string("Pull_") + hist->GetName()).c_str(), "Pull", bins, x0, x1 );
          for (int i=0;i<hist->GetEntries();i++){
            double p=0;
            double d = hist->GetBinContent(i) - data_plot->GetBinContent(i);
            double s2 = hist->GetBinContent(i) + data_plot->GetBinContent(i);
            double x = hist->GetBinCenter(i);
            if (s2!=0) p=d/sqrt(s2);
            pull->Fill(x, p);
          }
          pull->Write();


        }
        auto p2 = sigMCevents_tag.eventType().defaultProjections(nBins);
        for( unsigned i = 0 ; i != p2.size() -1; ++i )
        {
          for( unsigned j=i+1; j < p2.size(); ++j )
          {
            auto dalitz = Projection2D( p2[i], p2[j] );
            auto hdalitz = dalitz.plot(prefix);
            auto data_plot = sigevents_tag.makeProjection(dalitz, Prefix(prefix));
            real_t integral = 0;
            for( unsigned event = 0 ; event != sigMCevents_tag.size(); ++event )
            {
              auto pos = dalitz(sigMCevents_tag[event]);
              //hdalitz->Fill( pos.first, pos.second, cs_tag.prob( sigMCevents_tag[event], tagMCevents_tag[event] ) );
              hdalitz->Fill( pos.first, pos.second, std::norm(cs_tag.getVal( sigMCevents_tag[event], tagMCevents_tag[event] ))/norm );
              integral += std::norm(cs_tag.getVal( sigMCevents_tag[event], tagMCevents_tag[event] ) )/norm;
            }
          auto data_integral = data_plot->Integral();
          //auto hdalitz_integral = hdalitz->Integral();
          auto hdalitz_integral = integral;
          INFO("data integral = "<<data_integral);
          INFO("hdalitz integral = "<<hdalitz_integral);
          INFO("hdalitz integral = "<<integral);

            //hdalitz->Scale( data_plot->Integral() / hdalitz->Integral() );
            hdalitz->Scale( data_plot->Integral() / integral );
            hdalitz->SetName( ( std::string("MC_") + hdalitz->GetName() ).c_str() );
            hdalitz->Write();
            data_plot->Write(); 

            auto px = hdalitz->ProjectionX();
            auto py = hdalitz->ProjectionY();

            auto bins = px->GetBin(px->GetEntries() -1);
            auto x0 = px->GetBinCenter(0);
            auto x1 = px->GetBinCenter(bins +1);

            bins = py->GetBin(py->GetEntries() -1);
            auto y0 = py->GetBinCenter(0);
            auto y1 = py->GetBinCenter(bins +1);


//            bins = int(sqrt(bins));


            TH2D * pull_2D = new TH2D( (std::string("Pull_") + hdalitz->GetName() ).c_str(), "Pull", bins, x0, x1, bins, y0, y1 );
            for (int i=0;i<data_plot->GetEntries();i++){
              for (int j=0;j<hdalitz->GetEntries();j++){
                double p=0;
                double d = hdalitz->GetBinContent(i,j) - data_plot->GetBinContent(i,j);
                double s2 = hdalitz->GetBinContent(i,j) + data_plot->GetBinContent(i,j);
                if (s2!=0) p = d/sqrt(s2);
                double x = px->GetBinCenter(i);
                double y = py->GetBinCenter(j);
                pull_2D->Fill(x, y, p);
            }
          }
          pull_2D->Write();
        }
        }
        f->Close();
        
}



        void setEvents(int i, EventList list1, EventList list2){
          m_Psi[i].setEvents(list1, list2);
        }
        void setMC(int i, EventList list1, EventList list2){
          m_Psi[i].setMC(list1, list2);
        }

        pCorrelatedSum pdf(int i){
          return m_Psi[i];
        }
        

        CombCorrLL(std::vector<EventList> SigData, 
                   std::vector<EventList> TagData, 
               
               
                   EventType SigType,
                   std::vector<EventType> TagType,
                   MinuitParameterSet mps,
                   size_t seed,
                   size_t NInt
        ):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        //m_SigInt(SigInt),
               
                        //m_SigType(SigType),
                        m_TagType(TagType),
                        m_mps(mps),
                       // m_SumFactors(sumFactors),
                        m_debug(NamedParameter<bool>("CombCorrLL::Debug", false, "Debug CombCorrLL"))
                        {
                              TRandom3 rndm;
                              rndm.SetSeed( seed );
                              gRandom = &rndm;


                           
                           

                            m_MC= Generator<>(SigType, &rndm).generate(NInt);
                            INFO("n(m_MC) = "<<m_MC.size());

                            std::vector<pCorrelatedSum> pCS = {};
                            for (auto i=0; i < m_SigData.size() ; i++){

                               m_TagInt.push_back(Generator<>(TagType[i], &rndm).generate(NInt));
                              if (m_debug) INFO("Making pCorrelatedSum");
                                pCorrelatedSum _pCS = pCorrelatedSum(SigType, TagType[i], mps);


                                //m_TagInt.push_back(mcTag);
                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                _pCS.setMC(m_MC, m_TagInt[i]);
                                INFO("Preparing "<<i+1<<"/"<<m_SigData.size());
                                _pCS.prepare();
                                INFO("Prepared "<<i+1<<"/"<<m_SigData.size());

                             
                              _pCS.debugNorm();
                                INFO("norm = "<<_pCS.norm());
                                m_Psi.push_back(_pCS);
                                INFO("Finished making LL"<<i);
                            }
                            INFO("Made LLs");
                            INFO("norm[0] = "<< m_Psi[0].norm());
    //                        m_Psi(*pCS);


                        }
     

        CombCorrLL(std::vector<EventList> SigData, 
                   std::vector<EventList> TagData, 
                   std::vector<EventList> SigInt, 
                   std::vector<EventList> TagInt,
                   std::vector<EventType> SigType,
                   std::vector<EventType> TagType,
                   MinuitParameterSet mps,
                   std::vector<std::string> sumFactors):
                        m_SigData(SigData),
                        m_TagData(TagData),
                        m_SigInt(SigInt),
                        m_TagInt(TagInt),
                        m_SigType(SigType),
                        m_TagType(TagType),
                        m_mps(mps),
                        m_SumFactors(sumFactors),
                        m_debug(NamedParameter<bool>("CombCorrLL::Debug", false, "Debug CombCorrLL"))
                        {
                            std::vector<pCorrelatedSum> pCS = {};
                            m_MC = SigInt[0];
                            for (auto i=0; i < m_SigData.size() ; i++){
                                pCorrelatedSum _pCS = pCorrelatedSum(m_SigType[0], m_TagType[i], m_mps, m_SumFactors[i]);
                                _pCS.setEvents(m_SigData[i], m_TagData[i]);
                                //_pCS.setMC(m_SigInt[i], m_TagInt[i]);
                                _pCS.setMC(m_SigInt[0], m_TagInt[i]);
                                if (m_debug) INFO("Preparing "<<i+1<<"/"<<m_SigData.size());

                                _pCS.prepare();
                                if (m_debug) INFO("Prepared "<<i+1<<"/"<<m_SigData.size());

                             // INFO("psiMC(0,0) = "<<_pCS.getVal(m_MC[0], m_TagInt[i][0]));
                              _pCS.debugNorm();
                                if (m_debug) INFO("Norm = "<<_pCS.norm());
                                if (m_debug) INFO("LL = "<<_pCS.LL());

                                m_Psi.push_back(_pCS);
                            }

    //                        m_Psi(*pCS);


                        }
        void reset(){
            for (size_t i=0;i<m_Psi.size();i++) m_Psi[i].reset(true);
        }


        double LL(int i){

                if (m_debug)INFO("At "<<i<<" out of "<<m_Psi.size());
                //auto _LL =  make_likelihood( m_SigData[i], m_TagData[i] , m_Psi[i]);
                auto psi = m_Psi[i];
                if (m_debug) INFO("Getting norm for "<<i);
              
               
              
                real_t norm = psi.norm();
                if (m_debug) INFO("LL = "<<psi.LL());
                if (m_debug) INFO("Norm = "<<norm);

                real_t ll = 0;
                #pragma omp parallel for reduction( +: ll )
                for ( unsigned int j = 0; j < m_SigData[i].size(); ++j ) {
                  auto evt1 = m_SigData[i][j];
                  auto evt2 = m_TagData[i][j];
                  auto prob = std::norm(psi.getVal(evt1,evt2))/norm;
                  ll += -2* log(prob);

                }
            if (m_debug) INFO("LL = "<<ll);

            return ll;

        }

        double getVal(){
          real_t ll=0;
          for (int i=0;i<m_Psi.size();i++) ll+=LL(i);

          return ll;

        }

  };
}
#endif
