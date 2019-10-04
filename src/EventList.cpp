#include <Rtypes.h>
#include <TDirectory.h>
#include <TEventList.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>

#include <stddef.h>
#include <complex>
#include <functional>
#include <iterator>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Projection.h"
#include "AmpGen/TreeReader.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Event.h"
#include "AmpGen/Types.h"
#include "AmpGen/ProfileClock.h"
using namespace AmpGen;

EventList::EventList( const EventType& type ) : m_eventType( type ) {} 

void EventList::loadFromFile( const std::string& fname, const ArgumentPack& args )
{
  auto current_file = gFile; 
  auto tokens = split( fname, ':'); 
  TTree* tree = nullptr;
  if( fname == "" ) FATAL("Filename must be specified to load data"); 
  if( tokens.size() == 2 ){
    gFile = TFile::Open( tokens[0].c_str(), "READ"); 
    if( gFile == nullptr ) FATAL("Failed to load file: " << tokens[0] );
    tree = (TTree*)gFile->Get( tokens[1].c_str() );
  }
  else {
    gFile = TFile::Open( fname.c_str(), "READ");
    if( gFile == nullptr ) FATAL("Failed to load file: " << tokens[0] );
    tree = (TTree*)gFile->Get("DalitzEventList");
  }
  if( tree == nullptr ) FATAL( "Failed to load tree from file: " << fname );
  loadFromTree( tree, args );
  gFile->Close();
  gFile = current_file; 
}

void EventList::loadFromTree( TTree* tree, const ArgumentPack& args )
{
  ProfileClock read_time; 
  auto pdfSize      = args.getArg<CacheSize>(0).val;
  auto filter       = args.getArg<Filter>(std::string("")).val;
  auto getGenPdf    = args.getArg<GetGenPdf>(false).val;
  auto weightBranch = args.getArg<WeightBranch>(std::string("")).val;
  auto branches     = args.getArg<Branches>().val;
  auto applySym     = args.getArg<ApplySym>(false).val;
  auto entryList    = args.getArg<EntryList>().val; 
  auto eventFormat  = m_eventType.getEventFormat( true );

  Event temp( branches.size() == 0 ? eventFormat.size() : branches.size() , pdfSize );
  temp.setWeight( 1 );
  temp.setGenPdf( 1 );
  tree->SetBranchStatus( "*", 0 );

  TreeReader<real_t> tr( tree );
  if( branches.size() != 0 ){
    INFO("Branches = [" << vectorToString(branches, ", ") << "]" );
    for ( auto branch = branches.begin(); branch != branches.end(); ++branch ) {
      unsigned int pos = std::distance( branches.begin(), branch );
      tr.setBranch( *branch, &(temp[pos]) );
      if( pos >= eventFormat.size() ){
        INFO("Specifiying event extension: " << *branch << " " << pos << " " << eventFormat.size() );
        m_extensions[ *branch ] = pos; 
      }
    }
  }
  else {
    for ( auto& branch : eventFormat ){
      tr.setBranch( branch.first, &(temp[branch.second]) );
    }
  }
  if( getGenPdf )          tr.setBranch( "genPdf",     temp.pGenPdf() );
  if( weightBranch != "" ) tr.setBranch( weightBranch, temp.pWeight() );
  if( filter != "" ){
    if( entryList.size() != 0 ){
      WARNING("Specified entry list and filter, will overwrite list with specified selection");
    }
    tr.prepare();
    tree->Draw(">>evtList", filter.c_str() );
    TEventList* evtList = (TEventList*)gDirectory->Get("evtList");
    for( int i = 0 ; i < evtList->GetN(); ++i ) 
      entryList.push_back( evtList->GetEntry(i) );
  }
  bool hasEventList    = entryList.size() != 0;
  unsigned int nEvents = hasEventList ? entryList.size() : tree->GetEntries();
  m_data.reserve( nEvents );
  auto symmetriser = m_eventType.symmetriser();
  for ( unsigned int evt = 0; evt < nEvents; ++evt ) {
    tr.getEntry( hasEventList ? entryList[evt] : evt );
    if( applySym ) symmetriser( temp );
    m_data.push_back( temp );
  }
  read_time.stop();
  INFO("Time to read tree = " << read_time << "[ms]; nEntries = " << size() );
}

TTree* EventList::tree( const std::string& name, const std::vector<std::string>& extraBranches )
{
  TTree* outputTree = new TTree( name.c_str(), name.c_str() );
  if ( size() == 0 ) {
    ERROR( "Trying to output empty tree" );
    return nullptr;
  }
  Event tmp = *( begin() );
  double genPdf = 1;
  double weight = 1;
  auto format = m_eventType.getEventFormat( true );
  
  for ( auto& f : format ){
    outputTree->Branch( f.first.c_str(), tmp.address( f.second ) );
  }
  for ( auto& f : m_extensions ){
    outputTree->Branch( f.first.c_str(), tmp.address( f.second ) );
  } 
  outputTree->Branch( "genPdf", &genPdf );
  outputTree->Branch( "weight", &weight );
  for ( auto& evt : *this ) {
    tmp    = evt;
    genPdf = evt.genPdf();
    weight = evt.weight();
    outputTree->Fill();
  }
  return outputTree;
}

std::vector<TH1D*> EventList::makeProjections( const std::vector<Projection>& projections, const ArgumentPack& args )
{
  std::vector<TH1D*> plots;
  for ( auto& proj : projections ) {
    TH1D* plot = makeProjection(proj, args );
    DEBUG("Made plot ... " << plot->GetName() );
    plots.push_back( plot );
  }
  return plots;
}

TH1D* EventList::makeProjection( const Projection& projection, const ArgumentPack& args ) const 
{
  auto selection      = args.getArg<Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<Prefix>(std::string(""));
  auto plot = projection.plot(prefix);
  plot->SetLineColor(args.getArg<LineColor>(kBlack).val); 
  plot->SetMarkerSize(0);
  for( auto& evt : m_data ){
    if( selection != nullptr && !selection(evt) ) continue;
    auto pos = projection(evt);
    plot->Fill( pos, evt.weight() * ( weightFunction == nullptr ? 1 : weightFunction(evt) / evt.genPdf() ) );
  }
  if( selection != nullptr ) INFO("Filter efficiency = " << plot->GetEntries() << " / " << size() );
  return plot;
}

TH2D* EventList::makeProjection( const Projection2D& projection, const ArgumentPack& args ) const
{
  auto selection      = args.getArg<Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<Prefix>().val;
  auto plot           = projection.plot(prefix);
  for ( auto& evt : m_data ){
    if ( selection != nullptr && !selection(evt) ) continue;
    auto pos = projection(evt);
    plot->Fill( pos.first, pos.second, evt.weight() * ( weightFunction == nullptr ? 1 : weightFunction(evt) / evt.genPdf() ) );
  }
  return plot;
}

void EventList::printCacheInfo( const unsigned int& nEvt )
{
  for ( auto& ind : m_pdfIndex ) {
    INFO( "Cache[" << ind.second << "] = " << ind.first << " = " << at( nEvt ).getCache( ind.second ) );
  }
}

size_t EventList::getCacheIndex( const CompiledExpressionBase& PDF ) const
{
  auto pdfIndex = m_pdfIndex.find( FNV1a_hash( PDF.name() ) );
  if ( pdfIndex != m_pdfIndex.end() )
    return pdfIndex->second;
  else
    ERROR( "FATAL: PDF Index for " << PDF.name() << " not found" );
  return 999;
}

size_t EventList::getCacheIndex( const CompiledExpressionBase& PDF, bool& isRegistered ) const
{
  auto pdfIndex = m_pdfIndex.find( FNV1a_hash( PDF.name() ) );
  if ( pdfIndex != m_pdfIndex.end() ) {
    isRegistered = true;
    return pdfIndex->second;
  }
  isRegistered = false;
  return 999;
}

void EventList::resetCache()
{
  m_pdfIndex.clear();
  for ( auto evt = begin(); evt != end(); ++evt ) evt->resizeCache( 0 );
  m_lastCachePosition = 0;
}

double EventList::integral() const
{
  double integral = 0;
  for ( auto& evt : *this ) {
    integral += evt.weight();
  }
  return integral;
}
void EventList::add( const EventList& evts )
{
  resetCache();
  WARNING( "Adding event lists invalidates cache state" );
  for ( auto& evt : evts ) {
    m_data.push_back( evt );
    rbegin()->resizeCache( 0 );
  }
}
double EventList::norm()
{
  if ( m_norm == 0 ) {
    double totalWeight = 0;
#pragma omp parallel for reduction( + : totalWeight )
    for ( unsigned int i = 0; i < size(); ++i ) totalWeight += ( *this )[i].weight() / ( *this )[i].genPdf();
    m_norm               = totalWeight;
  }
  return m_norm;
}

void EventList::clear() 
{ 
  m_data.clear(); 
}

void EventList::erase(const std::vector<Event>::iterator& begin, 
                      const std::vector<Event>::iterator& end)
{
  m_data.erase( begin, end );
}

void EventList::reserveCache(const size_t& size)
{ 
  if ( size >= at(0).cacheSize() )
    for (auto& evt : *this) evt.resizeCache(evt.cacheSize() + size);
}

void EventList::resizeCache(const size_t& newCacheSize )
{
  for (auto& evt : *this) evt.resizeCache( newCacheSize );
}

