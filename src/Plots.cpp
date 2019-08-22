#if __cplusplus >= 201402L
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Plots.h"
#include "AmpGen/DalitzIntegrator.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <chrono>
using namespace AmpGen;
void AmpGen::perAmplitudePlot( const EventList& evts, 
    const Projection& projection,
    const CoherentSum& pdf )
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
      auto pdf_i             = pdf[i].amp;
      auto pdf_j             = pdf[j].amp;
      unsigned int index_i   = evts.getCacheIndex(pdf_i);
      unsigned int index_j   = evts.getCacheIndex(pdf_j);
      const std::string name = pdf_i.name() + "_" + pdf_j.name();
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

TGraph* AmpGen::boundary(const AmpGen::EventType& eventType, 
                         const std::function<double(const AmpGen::Event&)>& p1, 
                         const std::function<double(const AmpGen::Event&)>& p2 )
{
  auto s0 = pow(eventType.motherMass(),2); 
  auto s1 = pow(eventType.mass(0),2); 
  auto s2 = pow(eventType.mass(1),2); 
  auto s3 = pow(eventType.mass(2),2); 

  DalitzIntegrator di( s0, s1, s2, s3 );

  TGraph* gboundary = new TGraph();

  Event tmp(12);

  for( double x = 0 ; x <= 1; x+=0.001){
    di.setEvent( {x,0}, tmp.address()  );
    gboundary->SetPoint( gboundary->GetN(), p1(tmp), p2(tmp) ); 
  }
  for( double y = 0 ; y <= 1; y+=0.01){
    di.setEvent( {1,y}, tmp.address()  );
    gboundary->SetPoint( gboundary->GetN(), p1(tmp), p2(tmp) ); 
  }
  for( double x = 0 ; x <= 1; x+=0.001){
    di.setEvent( {1-x,1}, tmp.address()  );
    gboundary->SetPoint( gboundary->GetN(), p1(tmp), p2(tmp) ); 
  }
  for( double y = 0 ; y <= 1; y+=0.01){
    di.setEvent( {0,1-y}, tmp.address()  );
    gboundary->SetPoint( gboundary->GetN(), p1(tmp), p2(tmp) ); 
  }
  return gboundary;
} 

#endif
