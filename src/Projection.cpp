#include "AmpGen/Projection.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include <stdio.h>

#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"

using namespace AmpGen;
using namespace AmpGen::PlotOptions; 

Projection::Projection() = default; 

Projection::Projection( const std::function<double(const Event&)>& fcn, 
    const std::string& name, const std::string& xAxisTitle,
    const size_t& nBins, const double& min, const double& max,
    const std::string& units 
    ) :
  m_func(fcn),
  m_name(name),
  m_xAxisTitle(xAxisTitle),
  m_units(units),
  m_nBins(nBins), 
  m_min(min),
  m_max(max)
{
  m_width = (m_max-m_min)/double(m_nBins);
}

TH1D* Projection::plot(const std::string& prefix) const {
  std::string p = ( prefix==""?"":prefix+"_");
  TH1D* plot = new TH1D( (p + m_name ).c_str(),"",m_nBins,m_min,m_max);
  if( m_units != "" ){ 
    plot->GetXaxis()->SetTitle( mysprintf("%s \\left[%s\\right]", m_xAxisTitle.c_str(),m_units.c_str() ).c_str() );
    plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f %s)", m_width,m_units.c_str()).c_str() );
  }
  else {
    plot->GetYaxis()->SetTitle( mysprintf("\\mathrm{Entries} / (%0.2f)",m_width).c_str() );
    plot->GetXaxis()->SetTitle( m_xAxisTitle.c_str() );
  }
  plot->GetYaxis()->SetTitleOffset(1.35);
  plot->SetMarkerSize(0);
  plot->SetMinimum(0);
  
  DEBUG("Returning plot: [" << m_min << " " << m_max << "] " << m_name << " " << 
     plot->GetXaxis()->GetBinLowEdge(1) << " " <<
     plot->GetXaxis()->GetBinLowEdge(1 + m_nBins)
      );
  return plot;
}
std::function<int(const Event& evt )> Projection::binFunctor() const {
  return [this](auto& evt){ return int ( ( (*this)(evt) - m_min ) / m_width ) ;};
}

TH2D* Projection2D::plot(const std::string& prefix) const {
  std::string p = ( prefix==""?"":prefix+"_");
  TH2D* plot = new TH2D( ( p + xAxis.m_name +"_vs_"+yAxis.m_name).c_str(),"",
      xAxis.m_nBins,xAxis.m_min,xAxis.m_max ,
      yAxis.m_nBins,yAxis.m_min,yAxis.m_max );

  plot->GetXaxis()->SetTitle( xAxis.m_xAxisTitle.c_str() ); 
  plot->GetYaxis()->SetTitle( yAxis.m_xAxisTitle.c_str() ); 
  plot->GetYaxis()->SetTitleOffset(1.35);
  return plot; 
}

const std::string Projection::name() const { return m_name; }
double Projection::operator()( const Event& evt ) const { return m_func( evt ); }

std::pair<double, double> Projection2D::operator()( const Event& evt ) const
{
  return {xAxis.m_func( evt ), yAxis.m_func( evt )};
}

template <> TH1D* Projection::projInternal( const EventList& events, const ArgumentPack& args) const 
{ 
  auto selection      = args.getArg<PlotOptions::Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<PlotOptions::Prefix>(std::string(""));
  auto axis           = plot(prefix);
  axis->SetLineColor(args.getArg<PlotOptions::LineColor>(kBlack).val); 
  axis->SetMarkerSize(0);
  for( auto& evt : events )
  {
    if( selection != nullptr && !selection(evt) ) continue;
    auto pos = operator()(evt);
    axis->Fill( pos, evt.weight() * ( weightFunction == nullptr ? 1 : weightFunction(evt) / evt.genPdf() ) );
  }
  if( selection != nullptr ) INFO("Filter efficiency = " << axis->GetEntries() << " / " << events.size() );
  return axis;
}

template <> std::tuple<std::vector<TH1D*>, THStack*> Projection::projInternal(const EventList& events, const Projection::keyedFunctors& weightFunction, const ArgumentPack& args) const
{
  std::vector<TH1D*> hists; 
  double norm_sum = args.getArg<Norm>(1).val;
  std::string prefix = args.getArg<PlotOptions::Prefix>().val;
  bool autowrite     = args.get<PlotOptions::AutoWrite>() != nullptr;
  THStack* stack     = args.getArg<PlotOptions::AddTo>(new THStack()).val;
  auto selection      = args.getArg<Selection>().val;
  if( prefix != "" ) prefix = prefix +"_";
  for( auto& key : weightFunction.keys ) 
    hists.push_back( plot(prefix + key ) );
  for( const auto& evt : events ){
    if( selection != nullptr && !selection(evt) ) continue;
    auto pos = operator()(evt);
    auto weights = weightFunction(evt);
    for( unsigned j = 0 ; j != weights.size(); ++j ) hists[j]->Fill( pos, evt.weight() * weights[j] / evt.genPdf() ); 
  }
  std::sort( std::begin(hists), std::end(hists), [](auto& h1, auto& h2){ return h1->Integral() < h2->Integral() ; } );
  double total = std::accumulate( std::begin(hists), std::end(hists), 0.0, [](double& t, auto& h){ return t + h->Integral() ; } ); 
  if( total == 0 ) ERROR("Norm = " << total );
  else for( auto& h : hists ) h->Scale( norm_sum / total );
  stack->SetName( (prefix + name() + "_stack").c_str());
  for( auto& h : hists ){
    stack->Add(h, "C HIST");
    if( autowrite ) h->Write();
  }
  if( autowrite ) stack->Write();
  return {hists, stack};
}

#if ENABLE_AVX
template <> TH1D* Projection::projInternal( const EventListSIMD& events, const ArgumentPack& args) const 
{ 
  return events.makeProjection(*this, args); 
}

template <> std::tuple<std::vector<TH1D*>, THStack*> Projection::projInternal(const EventListSIMD& events, const Projection::keyedFunctors& weightFunction, const ArgumentPack& args) const
{
  std::vector<TH1D*> hists; 
  double norm_sum = args.getArg<Norm>(1).val;
  std::string prefix = args.getArg<PlotOptions::Prefix>().val;
  bool autowrite     = args.get<PlotOptions::AutoWrite>() != nullptr;
  THStack* stack     = args.getArg<PlotOptions::AddTo>(new THStack()).val;
  auto selection      = args.getArg<Selection>().val;
  if( prefix != "" ) prefix = prefix +"_";
  for( auto& key : weightFunction.keys ) 
    hists.push_back( plot(prefix + key ) );
  for( const auto& evt : events ){
    if( selection != nullptr && !selection(evt) ) continue;
    auto pos = operator()(evt);
    auto weights = weightFunction(evt);
    for( unsigned j = 0 ; j != weights.size(); ++j ) hists[j]->Fill( pos, evt.weight() * weights[j] / evt.genPdf() ); 
  }
  std::sort( std::begin(hists), std::end(hists), [](auto& h1, auto& h2){ return h1->Integral() < h2->Integral() ; } );
  double total = std::accumulate( std::begin(hists), std::end(hists), 0.0, [](double& t, auto& h){ return t + h->Integral() ; } ); 
  if( total == 0 ) ERROR("Norm = " << total );
  else for( auto& h : hists ) h->Scale( norm_sum / total );
  stack->SetName( (prefix + name() + "_stack").c_str());
  for( auto& h : hists ){
    stack->Add(h, "C HIST");
    if( autowrite ) h->Write();
  }
  if( autowrite ) stack->Write();
  return {hists, stack};
}
#endif
