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

template <typename key_type, typename value_type, typename default_type=key_type> key_type inv_map(const std::map<key_type, value_type>& mp, const value_type& k, const default_type& df = default_type()) 
{
  for( const auto& [key,val] : mp ) if( val == k ) return key; 
  return df;
}

void EventList::loadFromTree( TTree* tree, const ArgumentPack& args )
{
  ProfileClock read_time; 
  if( m_eventType.size() == 0 ){
    auto tokens = split( tree->GetTitle(), ' ');
    if( tokens.size() != 1 ) m_eventType = EventType( tokens ); 
    INFO("Attempted automatic deduction of eventType: " << m_eventType );
  } 
  auto filter       = args.getArg<Filter>(std::string("")).val;
  auto getGenPdf    = args.getArg<GetGenPdf>(false).val;
  auto weightBranch = args.getArg<WeightBranch>(std::string("")).val;
  auto branches     = args.getArg<Branches>().val;
  auto extraBranches= args.getArg<ExtraBranches>().val; 
  auto applySym     = args.getArg<ApplySym>(false).val;
  auto entryList    = args.getArg<EntryList>().val; 
  auto eventFormat  = m_eventType.getEventFormat( true );
  auto inputUnits   = args.getArg<InputUnits>(Units::GeV);
  auto idBranches   = args.getArg<IdBranches>({}).val;
  Event temp( eventFormat.size() + extraBranches.size());
  temp.setWeight( 1 );
  temp.setGenPdf( 1 );
  tree->SetBranchStatus( "*", 0 );
  TreeReader tr( tree );
  bool hasEnergy = branches.size() == 0 || branches.size() == 4 * m_eventType.size(); // if the energy of the particle has been explicitly specified //   
  std::vector<int> ids( m_eventType.size()  ); 
  if( branches.size() != 0 )
  {
    DEBUG("Branches = [" << vectorToString(branches, ", ") << "]" );
    for (unsigned p = 0 ; p != branches.size(); ++p ) 
    {
      auto pos = hasEnergy ? p : 4 * int(p/3) + p % 3 ;  
      DEBUG("Setting branch: " << branches[p] << " pos: " << pos << " fmt = " << inv_map( eventFormat, pos, "NOT FOUND" ) << " has energy? " << hasEnergy );
      tr.setBranch( branches[p], &(temp[pos]) );
    }
    if( idBranches.size() != 0  )
    {
      if( idBranches.size() != m_eventType.size() ) FATAL("Number of ID branches should be number of final state particles"); 
      for( int i = 0; i != ids.size(); ++i ) tr.setBranch( idBranches[i], ids.data() + i);
    }
  }
  else for ( auto& branch : eventFormat ) tr.setBranch( branch.first, &(temp[branch.second]) );
  auto pos = eventFormat.size();
  for( const auto& branch : extraBranches ){ 
    tr.setBranch( branch, &(temp[pos]) );
    m_extensions[branch] = pos++;
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
  if( entryList.size() != 0 ) tr.setEntryList( entryList );
  m_data.reserve( tr.nEntries() );
  auto symmetriser = m_eventType.symmetriser();
  auto automaticOrdering = m_eventType.automaticOrdering();
  for (const auto& evt : tr) {
    if( inputUnits != Units::GeV ) for( unsigned k = 0; k != eventFormat.size(); ++k ) temp[k] *= to_double(inputUnits); 
    if( idBranches.size() != 0 && !automaticOrdering(temp, ids) ) 
      WARNING("Failed to order event: " << evt ); 
    if( applySym ) symmetriser(temp); 
    if( ! hasEnergy ){ 
      for( unsigned int k = 0 ; k != m_eventType.size(); ++k ) 
        temp[4*k + 3] = sqrt( m_eventType.mass(k) * m_eventType.mass(k) + temp[4*k+0]*temp[4*k+0] + temp[4*k+1]*temp[4*k+1] + temp[4*k+2]*temp[4*k+2] ); 
    }
    push_back( temp );
  }
  read_time.stop();
  INFO("Time to read tree = " << read_time << "[ms]; nEntries = " << size() );
}

TTree* EventList::tree( const std::string& name, const std::vector<std::string>& extraBranches ) const
{
  std::string title = m_eventType.mother();
  for( unsigned i = 0 ; i != m_eventType.size(); ++i ) title += " " + m_eventType[i]; 
  TTree* outputTree = new TTree( name.c_str(), title.c_str() );
  if ( size() == 0 ) {
    ERROR( "Trying to output empty tree" );
    return nullptr;
  }
  Event tmp = *( begin() );
  double genPdf = 1;
  double weight = 1;
  auto format = m_eventType.getEventFormat( true );
  for ( const auto& f : format ) outputTree->Branch( f.first.c_str(), tmp.address( f.second ) );  
  for ( const auto& f : m_extensions ) outputTree->Branch( f.first.c_str(), tmp.address( f.second ) );
  outputTree->Branch( "genPdf", &genPdf );
  outputTree->Branch( "weight", &weight );
  for ( const auto& evt : *this ) {
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
  for ( const auto& proj : projections ) plots.push_back( makeProjection(proj, args) );
  return plots;
}

TH1D* EventList::makeProjection( const Projection& projection, const ArgumentPack& args ) const 
{
  auto selection      = args.getArg<PlotOptions::Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<PlotOptions::Prefix>(std::string(""));
  auto plot = projection.plot(prefix);
  plot->SetLineColor(args.getArg<PlotOptions::LineColor>(kBlack).val); 
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
  auto selection      = args.getArg<PlotOptions::Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<PlotOptions::Prefix>().val;
  auto plot           = projection.plot(prefix);
  for ( auto& evt : m_data ){
    if ( selection != nullptr && !selection(evt) ) continue;
    auto pos = projection(evt);
    plot->Fill( pos.first, pos.second, evt.weight() * ( weightFunction == nullptr ? 1 : weightFunction(evt) / evt.genPdf() ) );
  }
  return plot;
}

double EventList::integral() const
{
  return std::accumulate( std::begin(*this), std::end(*this), 0, [](double rv, const auto& evt){ return rv + evt.weight(); } );
}

void EventList::add( const EventList& evts )
{
  for ( auto& evt : evts ) push_back( evt );
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

void EventList::reserve( const size_t& size ) 
{ 
  m_data.reserve( size ); 
}

void EventList::resize ( const size_t& size ) 
{ 
  m_data.resize(size); 
  for( unsigned int i = 0 ; i != size; ++i ) m_data[i].setIndex(i) ; 
} 

void EventList::push_back( const Event& evt ) 
{ 
  m_data.push_back( evt ); 
  m_data.rbegin()->setIndex(m_data.size()-1);
}

void EventList::emplace_back( const Event& evt) 
{ 
  m_data.emplace_back(evt) ; 
  m_data.rbegin()->setIndex(m_data.size()-1);
}
