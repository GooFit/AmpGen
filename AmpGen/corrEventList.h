/* Correlated Event - 2 events stuck together */
#ifndef CORREventList
#define CORREventList

#include "AmpGen/EventList.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
namespace AmpGen{

    class corrEventList{

        public:
            corrEventList(EventList& events1, EventList& events2);
            std::vector<TH1D *> makeProjection(const Projection& projection1, const Projection& projection2, const ArgumentPack& args) const;
            std::vector<std::vector<TH1D*> > makeProjections(const std::vector<Projection>& projections1, const std::vector<Projection>& projections2, const ArgumentPack& args) const;
            std::vector<TH2D *> makeProjection(const Projection2D& projection, const ArgumentPack& args) const;
            size_t size(){
                return m_events1.size();
            }
            EventType EventType1() const {return m_eventType1;}
            EventType EventType2() const {return m_eventType2;}
            template <class... ARGS> std::vector<std::vector<TH1D*> > makeDefaultProjections( const ARGS&... args )
            {
                auto argPack = ArgumentPack( args... );
                size_t nBins = argPack.getArg<Bins>(10000);
                std::cout<<"making projection with "<<nBins<<" bins\n";
                auto proj1 = m_eventType1.defaultProjections(nBins); 
                auto proj2 = m_eventType2.defaultProjections(nBins); 

                INFO("Projection1 has "<<proj1.size()<<"entries");
                INFO("Projection2 has "<<proj2.size()<<"entries");

                return makeProjections( proj1, proj2, argPack );
            }
            void setDebug(){
                m_debug=true;
            } 
        private:
            EventList& m_events1;
            EventList& m_events2;
            EventType m_eventType1;
            EventType m_eventType2;
            bool m_debug;


    };

  DECLARE_ARGUMENT(corrLineColor, int);
  DECLARE_ARGUMENT(corrDrawStyle, std::string);
  DECLARE_ARGUMENT(corrSelection, std::function<bool(const Event&, const Event& )>);
  DECLARE_ARGUMENT(corrWeightFunction, std::function<double(const Event&, const Event&)>);
  DECLARE_ARGUMENT(corrNorm, double);
  DECLARE_ARGUMENT(corrBranches, std::vector<std::string>);
  DECLARE_ARGUMENT(corrEntryList, std::vector<size_t>);
  DECLARE_ARGUMENT(corrGetGenPdf, bool);
  DECLARE_ARGUMENT(corrCacheSize, size_t);
  DECLARE_ARGUMENT(corrFilter, std::string);
  DECLARE_ARGUMENT(corrWeightBranch, std::string);      
  DECLARE_ARGUMENT(corrApplySym, bool);  
  DECLARE_ARGUMENT(corrPrefix, std::vector<std::string>);

} //namespace AmpGen





#endif