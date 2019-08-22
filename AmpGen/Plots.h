#ifndef AMPGEN_PLOTS_H
#define AMPGEN_PLOTS_H
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventList.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <chrono>

namespace AmpGen
{
  void perAmplitudePlot(const EventList& evts, const Projection& projection, const CoherentSum& pdf);

  template <size_t NBINS, size_t NROLLS>
    std::array<Bilinears, NBINS> getNorms( CoherentSum& fcn, BinnedIntegrator<NBINS, NROLLS>& bid )
    {

      std::array<Bilinears, NBINS> normalisations;
      for ( unsigned int i = 0; i < NBINS; ++i ) normalisations[i] = Bilinears( fcn.size(), fcn.size() );

      for ( unsigned int i = 0; i < fcn.size(); ++i ) {
        for ( unsigned int j = i; j < fcn.size(); ++j ) {
          bid.addIntegral( fcn[i].amp, fcn[j].amp, [i, j, &normalisations]( const auto& val ) {
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
        bid.addIntegral( fcn[i].amp, fcn[i].amp, [i, &normalisations]( const auto& val ) {
            for ( unsigned int bin = 0; bin < NBINS; ++bin ) normalisations[bin].set( i, 0, val[bin] );  
            } );
      }
      bid.flush();
      return normalisations;
    }

  template <size_t NBINS, class FCN>
    TH1D* plotWithError( EventList& events, FCN& fcn, const Projection& projection, const std::string& prefix,
        LinearErrorPropagator& linProp, const std::function<bool( const Event& )>& selection = nullptr )
    {
      BinnedIntegrator<NBINS, 10> bid( &events );
      if ( selection != nullptr ) bid.setSlice( selection );
      bid.setView( projection.binFunctor() );
      TH1D* plot = projection.plot();
      plot->SetName( ( prefix + plot->GetName() ).c_str() );
      auto normalisations = getNorms<NBINS>( fcn, bid );
      auto vectorBinFunctor = [&normalisations, &fcn, &bid] {
        fcn.transferParameters();
        bid.update( fcn, normalisations );
        std::vector<double> values(NBINS);
        double total = 0;
        for ( size_t bin = 0; bin < NBINS; ++bin ) {
          values[bin] = fcn.norm( normalisations[bin] );
          total += values[bin];
        }
        for ( size_t bin = 0; bin < NBINS; ++bin ) values[bin] /= total;
        return values;
      };
      auto values = vectorBinFunctor();
      auto errors = linProp.getVectorError( vectorBinFunctor, NBINS );
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
  TGraph* boundary(const AmpGen::EventType& type, 
                   const std::function<double(const AmpGen::Event&)>& p1, 
                   const std::function<double(const AmpGen::Event&)>& p2 );
} // namespace AmpGen

#endif
