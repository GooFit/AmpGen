#ifndef AMPGEN_CHI2CORRELATED_H
#define AMPGEN_CHI2CORRELATED_H

#include <functional>
#include <string>
#include <vector>

#include "AmpGen/BinDT.h"
#include <memory>
#include <ostream>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Event.h"

namespace AmpGen {





    class Chi2Correlated{
        struct Moment {
        double x;
        double xx;
        double N;
        std::vector<double> values;
        Moment() : x( 0 ), xx( 0 ), N( 0 ) {}
        void add( const double& value )
        {
            x += value;
            xx += value * value;
            N++;
            values.push_back( value );
        }
        void rescale( const double& val )
        {
            x *= val;
            xx *= ( val * val );
        }
        double val() { return x; }
        double var() { return N == 0 ? 0 : xx; }
        };
        public:
            Chi2Correlated( const EventList& dataEvents_sig, const EventList& dataEvents_tag,
                    const EventList& mcEvents_sig, const EventList& mcEvents_tag,
                    const std::function<double(const Event&, const Event&)>& fcn, 
                    const unsigned int& minEvents = 10):
                        m_binning_sig(dataEvents_sig, MinEvents(minEvents), Dim(dataEvents_sig.eventType().dof())),
                        m_binning_tag(dataEvents_tag, MinEvents(minEvents), Dim(dataEvents_tag.eventType().dof()))
                        {
                            doChi2(dataEvents_sig, dataEvents_tag, mcEvents_sig, mcEvents_tag, fcn);
                        }
            void doChi2(const EventList& dataEvents_sig, const EventList& dataEvents_tag,
                        const EventList& mcEvents_sig, const EventList& mcEvents_tag,
                        const std::function<double(const Event&, const Event&)>& fcn)
                        {
                            std::vector<Moment> data_sig(m_binning_sig.size());
                            std::vector<Moment> mc_sig(m_binning_sig.size());
                            std::vector<Moment> data_tag(m_binning_tag.size());
                            std::vector<Moment> mc_tag(m_binning_tag.size());
                            INFO("Splitting "<<dataEvents_sig.size()<<" data, "<<mcEvents_sig.size()<<" MC into "<<m_binning_sig.size()<<" bins");
                            unsigned int j=0;
                            double total_data_weight = 0;
                            double total_int_weight = 0;
                            for (int i=0;i<dataEvents_sig.size();i++){
                                auto d_sig = dataEvents_sig[i];
                                auto d_tag = dataEvents_tag[i];
                                if ( j % 1000000 == 0 && j != 0) INFO("Binned "<<j<<" data Events");
                                double w_sig = d_sig.weight();
                                double w_tag = d_tag.weight();
                                double w = w_sig * w_tag;
                                data_sig[m_binning_sig.getBinNumber(d_sig)].add(d_sig.weight());
                                data_tag[m_binning_tag.getBinNumber(d_tag)].add(d_tag.weight());
                                total_data_weight += w;

                                j++;
                            }
                            j =0 ;
                            for (int i=0;i<mcEvents_sig.size();i++){
                                auto evt_sig = mcEvents_sig[i];
                                auto evt_tag = mcEvents_tag[i];
                                if ( j % 1000000 == 0 && j != 0) INFO("Binned "<<j<<" sim. Events");
                                double w = fcn(evt_sig, evt_tag) * evt_sig.weight() * evt_tag.weight()/(evt_sig.genPdf() * evt_tag.genPdf());
                                mc_sig[m_binning_sig.getBinNumber(evt_sig)].add(w);
                                mc_tag[m_binning_tag.getBinNumber(evt_tag)].add(w);
                                total_int_weight += w;
                                j++;
                            }
                            double chi2 = 0;
                            for (unsigned int i=0;i<m_binning_sig.size();i++){
                                mc_sig[i].rescale(total_data_weight/total_int_weight);
                                double delta = data_sig[i].val() - mc_sig[i].val();
                                double sigma = data_sig[i].val() + mc_sig[i].val();
                                double tChi2 = delta * delta/sigma;
                                chi2 += tChi2;
                            }
                            m_chi2 = chi2;
                            m_nBins = m_binning_sig.size();
 
                        }
            double chi2() const {return m_chi2;}
            double nBins() const {return m_nBins;}
            void writeBinningToFile(const std::string& fileName_sig, const std::string& fileName_tag){
                m_binning_sig.serialize(fileName_sig);
                m_binning_tag.serialize(fileName_tag);
            }

                        
        private:
            double m_chi2;
            double m_nBins;
            BinDT m_binning_sig;
            BinDT m_binning_tag;

            
    };
}
#endif