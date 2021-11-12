#include "AmpGen/corrEventList.h"

namespace AmpGen{


    corrEventList::corrEventList(EventList& events1, EventList& events2):
        m_events1 (events1),
        m_events2 (events2),
        m_eventType1(events1.eventType()),
        m_eventType2(events2.eventType())
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

        //correlated part

        std::vector<TH1D *> hists;

        auto n = m_events1.size(); 
        int i=0;
        while (i < n){
            auto evt1 = m_events1[i];
            auto evt2 = m_events2[i];
            
            if( selection != nullptr && !selection(evt1, evt2) ) continue;
            auto pos1 = projection1(evt1);
            auto pos2 = projection2(evt2);
            plot1->Fill(pos1, evt1.weight() * (weightFunction == nullptr ? 1: weightFunction(evt1, evt2)/evt1.genPdf() ) );
            plot2->Fill(pos2, evt2.weight() * (weightFunction == nullptr ? 1: weightFunction(evt1, evt2)/evt2.genPdf() ) );

            if( selection != nullptr ) INFO("Filter efficiency = " << plot1->GetEntries() << " / " << m_events1.size() );
            if( selection != nullptr ) INFO("Filter efficiency = " << plot2->GetEntries() << " / " << m_events2.size() );

            hists = std::vector<TH1D *>({plot1, plot2});
            return hists;
        

        }
    return hists;


    }

    std::vector<std::vector<TH1D*> > corrEventList::makeProjections(const std::vector<Projection>& projections1, const std::vector<Projection>& projections2, const ArgumentPack& args) const {
        std::vector<std::vector<TH1D*> > plots;
        for (auto proj1 : projections1){
            for (auto proj2 : projections2){
                std::vector<TH1D*> plot = makeProjection(proj1, proj2, args);
                plots.push_back(plot);
            }
        }
        return plots;
    }



} //namespace AmpGen