#include "AmpGen/TreeReader.h"
using namespace AmpGen; 

TreeReader::TreeReader( TTree* tree ) : m_tree( tree ) {}

void TreeReader::setEntryList( const std::vector<size_t>& entryList ){ m_entryList = entryList ; }

void TreeReader::unsetEntryList(){ m_entryList.clear() ; }

void TreeReader::getEntry( const unsigned int& entry )
{
  if(!m_ready ) prepare(); 
  m_tree->GetEntry( m_entryList.size() == 0 ? entry : m_entryList[entry] );
  for ( auto& branch : m_branches ) branch->transfer();
}

void TreeReader::prepare()
{
  for ( auto& branch : m_branches ) {
    m_tree->SetBranchStatus( branch->name.c_str(), "1" );
    m_tree->SetBranchAddress( branch->name.c_str(), branch->address() );
  }
  m_ready = true;
}

size_t TreeReader::nEntries() const { return m_entryList.size() == 0 ? m_tree->GetEntries() : m_entryList.size(); }  

TreeReader::~TreeReader()
{
  for ( auto& branch : m_branches ) delete branch;
}

TreeReader::Iterator TreeReader::begin() { getEntry(0); return Iterator( 0, this ); }
TreeReader::Iterator TreeReader::end()   { return Iterator( nEntries(), this ); }
