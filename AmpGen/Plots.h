#ifndef AMPGEN_PLOTS_H
#define AMPGEN_PLOTS_H
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventList.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <chrono>

namespace AmpGen
{
  template <class PDF>
    void perAmplitudePlot( const EventList& evts, 
                           const Projection& projection,
                           PDF pdf )
    {
      struct PlotIJ {
        unsigned int i;
        unsigned int j;
        TH1D* hist;
        std::complex<double> amp;
      };

      TDirectory* dir = (TDirectory*)gFile->Get( ("perAmp_"+projection.name()).c_str() );
      if( dir == nullptr )
      {
        gFile->mkdir(  ("perAmp_"+ projection.name() ).c_str() );
        dir = (TDirectory*)gFile->Get( ("perAmp_"+projection.name()).c_str() );
      } 
      dir->cd();

      std::vector<std::pair<const Event*, double>> eventData;

      std::vector<PlotIJ> tmpPlots( pdf.size() * ( pdf.size() + 1 ) / 2 );

      unsigned int s = 0;
      for ( unsigned int i = 0; i < pdf.size(); ++i ) {

        for ( unsigned int j = i; j < pdf.size(); ++j ) {
          auto pdf_i             = pdf[i].pdf;
          auto pdf_j             = pdf[j].pdf;
          unsigned int index_i   = evts.getCacheIndex( pdf[i].pdf );
          unsigned int index_j   = evts.getCacheIndex( pdf[j].pdf );
          const std::string name = pdf_i.name() + "_" + pdf_j.name();
      //    INFO( name << " " << index_i << " " << index_j );
          tmpPlots[s].hist       = projection.plot(name);
          tmpPlots[s].i          = index_i;
          tmpPlots[s].j          = index_j;
          tmpPlots[s].amp        = pdf[i].coupling() * std::conj( pdf[j].coupling() );
          if ( index_i != index_j ) tmpPlots[s].amp = 2.0 * tmpPlots[s].amp;
          s++;
        }
      }
      for ( auto& evt : evts ) {
        double f = projection( evt );
        for ( auto& h : tmpPlots ) {
          std::complex<double> pdfValue = evt.getCache( h.i ) * std::conj( evt.getCache( h.j ) );
          double weight                 = std::real( h.amp * pdfValue ) * evt.weight() / evt.genPdf();
          h.hist->Fill( f, weight );
        }
      }
      for ( auto& h : tmpPlots ) {
        h.hist->Write();
        delete h.hist;
      }
      dir->Write();
      gFile->cd();
    }

  template <size_t NBINS, size_t NROLLS>
    std::array<Bilinears, NBINS> getNorms( CoherentSum& fcn, BinnedIntegrator<NBINS, NROLLS>& bid )
    {

      std::array<Bilinears, NBINS> normalisations;
      for ( unsigned int i = 0; i < NBINS; ++i ) normalisations[i] = Bilinears( fcn.size(), fcn.size() );

      for ( unsigned int i = 0; i < fcn.size(); ++i ) {
        for ( unsigned int j = i; j < fcn.size(); ++j ) {
          bid.addIntegral( fcn[i].pdf, fcn[j].pdf, [i, j, &normalisations]( const auto& val ) {
              for ( unsigned int bin = 0; bin < NBINS; ++bin ) {
              normalisations[bin].set( i, j, val[bin] );
              if ( i != j ) normalisations[bin].set( j, i, std::conj( val[bin] ) );
              }
              } );
        }
      }
      bid.flush();
      return normalisations;
    }

  template <size_t NBINS, size_t NROLLS>
    std::array<Bilinears, NBINS> getNorms( IncoherentSum& fcn, BinnedIntegrator<NBINS, NROLLS>& bid )
    {

      std::array<Bilinears, NBINS> normalisations;
      for ( unsigned int i = 0; i < NBINS; ++i ) normalisations[i] = Bilinears( fcn.size(), fcn.size() );

      for ( unsigned int i = 0; i < fcn.size(); ++i ) {
        bid.addIntegral( fcn[i].pdf, fcn[i].pdf, [i, &normalisations]( const auto& val ) {
            for ( unsigned int bin = 0; bin < NBINS; ++bin ) {
            normalisations[bin].set( i, 0, val[bin] );
            }
            } );
      }
      bid.flush();
      return normalisations;
    }

  template <size_t NBINS, class FCN>
    TH1D* plotWithError( EventList& events, FCN& fcn, const Projection& projection, const std::string& prefix,
        LinearErrorPropagator& linProp, const std::function<bool( const Event& )>& selection = nullptr )
    {

      bool hardcore = NamedParameter<bool>( "Hardcore", false );
      BinnedIntegrator<NBINS, 10> bid( &events );
      if ( selection != nullptr ) bid.setSlice( selection );
      bid.setView( projection.binFunctor() );

      TH1D* plot = projection.plot();
      plot->SetName( ( prefix + plot->GetName() ).c_str() );
      auto normalisations = getNorms<NBINS>( fcn, bid );

      auto vectorBinFunctor = [&normalisations, &fcn, &hardcore, &bid] {
        fcn.transferParameters();
        if ( hardcore ) bid.update( fcn, normalisations );
        std::array<double, NBINS> values;
        double total = 0;
        for ( size_t bin = 0; bin < NBINS; ++bin ) {
          values[bin] = fcn.norm( normalisations[bin] );
          total += values[bin];
        }
        for ( size_t bin = 0; bin < NBINS; ++bin ) values[bin] /= total;
        return values;
      };
      auto values = vectorBinFunctor();
      auto errors = linProp.getVectorError<NBINS>( vectorBinFunctor );
      for ( size_t bin = 0; bin < NBINS; ++bin ) {
        plot->SetBinContent( bin + 1, values[bin] );
        plot->SetBinError( bin + 1, errors[bin] );
      }
      return plot;
    }


  template <size_t NBINS, class FCN>
    std::vector<TH1D*> bandPlot( EventList& events, const std::string& prefix, FCN& fcn, LinearErrorPropagator& linProp )
    {
      std::vector<TH1D*> plots;
      auto axes = events.eventType().defaultProjections( NBINS );
      for ( auto& proj : axes ) {
        INFO( "Making plot:" << proj.name() );
        plots.push_back( plotWithError<NBINS>( events, fcn, proj, prefix, linProp ) );
      }
      return plots;
    }
} // namespace AmpGen

#endif
