#include "AmpGen/corrEventList.h"
#include "AmpGen/NamedParameter.h"
#include <fstream>
namespace AmpGen{


    corrEventList::corrEventList(EventList& events1, EventList& events2):
        m_events1 (events1),
        m_events2 (events2),
        m_eventType1(events1.eventType()),
        m_eventType2(events2.eventType()),
        m_debug(NamedParameter<bool>("corrEventList::debug", false, "Print Debug messages for corrEventList"))


    {




    }

    std::vector<TH1D *> corrEventList::makeProjection(const Projection& projection1, const Projection& projection2, const ArgumentPack& args) const{
        
                //From EventList
        auto selection      = args.getArg<corrSelection>().val;
        auto weightFunction = args.getArg<corrWeightFunction>().val;
        std::vector<std::string> prefix  = args.getArg<corrPrefix>(std::vector<std::string> {"", ""});
        auto plot1 = projection1.plot(prefix[0]);
        plot1->SetLineColor(args.getArg<LineColor>(kBlack).val); 
        plot1->SetMarkerSize(0);
        auto plot2 = projection2.plot(prefix[1]);
        plot2->SetLineColor(args.getArg<LineColor>(kBlack).val); 
        plot2->SetMarkerSize(0);

        //Make csv file with pos1,weight coordinates 
        std::ofstream file;
        std::stringstream stream;
        stream<<"test_"<<prefix[0]<<".csv";
        file.open(stream.str().c_str());

        //correlated part
        std::vector<TH1D *> pair;
        auto n = m_events1.size(); 
        if (m_debug) INFO("Projecting "<<n<<" pairs of events");
        int i=0;
        while (i < n){
            auto evt1 = m_events1[i];
            auto evt2 = m_events2[i];
            
            if( selection != nullptr && !selection(evt1, evt2) ) continue;
            auto pos1 = projection1(evt1);
            auto pos2 = projection2(evt2);
            if (i%1000==0){
            if (weightFunction!=nullptr) INFO("Weight Function = "<<weightFunction(evt1, evt2)/evt1.genPdf());
            INFO("Pos1 = "<<pos1*(-1));
            INFO("Pos2 = "<<pos2);
            }
            plot1->Fill(pos1);
            if (weightFunction!=nullptr){
                file<<pos1<<"\t"<<weightFunction(evt1, evt2)/evt1.genPdf()<<"\n";
            }
            //plot1->Fill(abs(pos1), weightFunction(evt1, evt2) );
            //plot1->Fill(abs(pos1), evt1.weight() * (weightFunction == nullptr ? 1: weightFunction(evt1, evt2)/evt1.genPdf() ) );
            plot2->Fill(pos2);
            
            //plot2->Fill(pos2,  weightFunction(evt1, evt2)/evt2.genPdf() );
            //plot2->Fill(pos2, evt2.weight() * (weightFunction == nullptr ? 1: weightFunction(evt1, evt2)/evt2.genPdf() ) );

            if( selection != nullptr ) INFO("Filter efficiency = " << plot1->GetEntries() << " / " << m_events1.size() );
            if( selection != nullptr ) INFO("Filter efficiency = " << plot2->GetEntries() << " / " << m_events2.size() );

            i++;

        }
            INFO("Mean of plot1 = "<<plot1->GetMean());
            pair.push_back(plot1);
            pair.push_back(plot2);
            file.close();

            return pair;

    }

    std::vector<std::vector<TH1D*> > corrEventList::makeProjections(const std::vector<Projection>& projections1, const std::vector<Projection>& projections2, const ArgumentPack& args) const {
        std::vector<std::vector<TH1D*> > plots;
        for (auto proj1 : projections1){
                std::vector<TH1D*> plot;
                std::cout<<"making projection1\n";
                if (projections2.size()<1){
                    plot = makeProjection(proj1, proj1, args);    
                }
                else{
                    for (auto proj2 : projections2){
                        std::cout<<"making projection2\n";
                    
                            plot = makeProjection(proj1, proj2, args);    
                    
                    }
                }
                plots.push_back(plot);
        }
        return plots;
    }



} //namespace AmpGen