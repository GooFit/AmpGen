#if ENABLE_AVX 

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
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Projection.h"
#include "AmpGen/TreeReader.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Event.h"
#include "AmpGen/Types.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/simd/utils.h"
using namespace AmpGen;

// ENABLE_DEBUG(EventListSIMD)

EventListSIMD::EventListSIMD( const EventType& type ) : 
  m_data(0, type.eventSize() ), 
  m_eventType( type ) {}

void EventListSIMD::loadFromFile( const std::string& fname, const ArgumentPack& args )
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
    if( tree == nullptr ) tree = (TTree*)gFile->Get("EventList");
  }
  if( tree == nullptr ) FATAL( "Failed to load tree from file: " << fname );
  loadFromTree( tree, args );
  gFile->Close();
  gFile = current_file; 
}

void EventListSIMD::loadFromTree( TTree* tree, const ArgumentPack& args )
{
  ProfileClock read_time; 
  if( m_eventType.size() == 0 ){
    auto tokens = split( tree->GetTitle(), ' ');
    if( tokens.size() != 1 ) setEventType( EventType( tokens ) ); 
    INFO("Attempted automatic deduction of eventType: " << m_eventType );
  } 
  auto filter       = args.getArg<Filter>(std::string("")).val;
  auto getGenPdf    = args.getArg<GetGenPdf>(false).val;
  auto weightBranch = args.getArg<WeightBranch>(std::string("")).val;
  auto branches     = args.getArg<Branches>().val;
  auto applySym     = args.getArg<ApplySym>(false).val;
  auto entryList    = args.getArg<EntryList>().val; 
  auto eventFormat  = m_eventType.getEventFormat( true );

  Event temp( branches.size() == 0 ? eventFormat.size() : branches.size());
  temp.setWeight( 1 );
  temp.setGenPdf( 1 );
  tree->SetBranchStatus( "*", 0 );

  TreeReader tr( tree );
  if( branches.size() != 0 ){
    INFO("Branches = [" << vectorToString(branches, ", ") << "]" );
    for ( auto branch = branches.begin(); branch != branches.end(); ++branch ) {
      unsigned int pos = std::distance( branches.begin(), branch );
      tr.setBranch( *branch, &(temp[pos]) );
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
    for( int i = 0 ; i < evtList->GetN(); ++i ) entryList.push_back( evtList->GetEntry(i) );
  }
  bool hasEventList    = entryList.size() != 0;
  size_t nEvents       = hasEventList ? entryList.size() : tree->GetEntries();
  std::array<Event, float_v::size> buffer;
  
  resize( nEvents ); 
  auto symmetriser = m_eventType.symmetriser();
  for ( unsigned int block = 0; block < m_data.nBlocks(); ++block ) 
  {
    for( unsigned k = 0 ; k != float_v::size; ++k )
    {
      auto evt = k + block * float_v::size; 
      if(evt  < m_data.size() )
      {
        tr.getEntry( hasEventList ? entryList[evt] : evt );
        if( applySym ) symmetriser( temp );
        buffer[k] = temp;  
      } 
      else buffer[k].setWeight(0);
    }
    gather( buffer, block );
  }
  read_time.stop();
  INFO("Time to read tree = " << read_time << "[ms]; nEntries = " << size() );
}


EventListSIMD::EventListSIMD( const EventList& other ) : EventListSIMD( other.eventType() ) 
{
  resize( other.size() );
  for( unsigned block = 0 ; block != m_data.nBlocks(); block++ )
  {
    for( unsigned j = 0 ; j != m_data.nFields(); ++j ) 
      m_data(block, j) = utils::gather<float_v>(other, [j](const auto& event){ return event[j]; } , block * float_v::size );
    m_weights[block] = utils::gather<float_v>(other,  [](const auto& event){ return event.weight(); }, block * float_v::size, 0);
    m_genPDF [block] = utils::gather<float_v>(other,  [](const auto& event){ return event.genPdf(); }, block * float_v::size, 1);
  }
} 



TTree* EventListSIMD::tree( const std::string& name, const std::vector<std::string>& extraBranches ) const
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
  // for ( const auto& f : m_extensions ) outputTree->Branch( f.first.c_str(), tmp.address( f.second ) );

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

std::vector<TH1D*> EventListSIMD::makeProjections( const std::vector<Projection>& projections, const ArgumentPack& args )
{
  std::vector<TH1D*> plots;
  for ( const auto& proj : projections ) plots.push_back( proj(*this,args) );
  return plots;
}

TH1D* EventListSIMD::makeProjection( const Projection& projection, const ArgumentPack& args ) const 
{
  return projection(*this, args);

}
TH2D* EventListSIMD::makeProjection( const Projection2D& projection, const ArgumentPack& args ) const
{
  auto selection      = args.getArg<PlotOptions::Selection>().val;
  auto weightFunction = args.getArg<WeightFunction>().val;
  std::string prefix  = args.getArg<PlotOptions::Prefix>().val;
  auto plot           = projection.plot(prefix);
  for ( const auto evt : *this ){
    if ( selection != nullptr && !selection(evt) ) continue;
    auto pos = projection(evt);
    plot->Fill( pos.first, pos.second, evt.weight() * ( weightFunction == nullptr ? 1 : weightFunction(evt) / evt.genPdf() ) );
  }
  return plot;
}


void EventListSIMD::clear() 
{ 
  m_data.clear(); 
}

const Event EventListSIMD::operator[]( const size_t& pos ) const 
{ 
  unsigned p = pos / float_v::size; 
  unsigned q = pos % float_v::size; 
  Event tempEvent( eventSize() );
  for( unsigned i = 0 ; i !=  tempEvent.size(); ++i ) 
    tempEvent[i] = m_data(p, i).at(q);
  tempEvent.setWeight( m_weights[p].at(q) );
  tempEvent.setGenPdf( m_genPDF[p].at(q) );
  tempEvent.setIndex( pos );
  return tempEvent; 
}

std::array<Event, AmpGen::float_v::size> EventListSIMD::scatter( unsigned pos ) const
{
  unsigned p = pos / float_v::size;
  std::array<Event, float_v::size> rt;
  auto vw = m_weights[p].to_array();
  auto vg = m_genPDF[p].to_array();
  for( unsigned evt = 0 ; evt != float_v::size; ++evt ){
    rt[evt] = Event( m_data.nFields() );
    rt[evt].setWeight(vw[evt]); 
    rt[evt].setGenPdf(vg[evt]); 
    rt[evt].setIndex(evt + pos);
  }
  for( unsigned field = 0 ; field != m_data.nFields(); ++field){
    auto v = m_data(p, field).to_array();
    for( unsigned evt = 0; evt != float_v::size; ++evt ) rt[evt][field] = v[evt]; 
  }
  return rt;
}

void EventListSIMD::gather( const std::array<Event, float_v::size>& data, unsigned pos )
{
  for( unsigned field = 0; field != m_data.nFields(); ++field ) 
    m_data(pos, field) = utils::gather<float_v>(data, [field](auto& event){ return event[field]; } );
  m_weights[pos] = utils::gather<float_v>(data, [](auto& event){ return event.weight() ; } );
  m_genPDF[pos]  = utils::gather<float_v>(data, [](auto& event){ return event.genPdf(); } );
}

#endif
