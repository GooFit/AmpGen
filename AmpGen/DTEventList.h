#include "AmpGen/DTEvent.h"
using namespace std;
using namespace AmpGen;
#ifndef DTEVENTLIST
#define DTEVENTLIST
struct DTEventList : public std::vector<DTEvent> 
{
  AmpGen::EventType m_sigType; 
  AmpGen::EventType m_tagType;
  DTEventList( const AmpGen::EventType& signal, const AmpGen::EventType& tag ) : m_sigType(signal), m_tagType(tag) {}
  std::string particleName(const AmpGen::EventType& type, const size_t& j);
  TTree* tree(const std::string& name);
};
#endif
std::string DTEventList::particleName(const AmpGen::EventType& type, const size_t& j)
{
  auto count = type.count(j);
  if( count.second == 1 ) return programatic_name(type[j]);
  return programatic_name(type[j])+std::to_string(count.first);
}

TTree* DTEventList::tree(const std::string& name)
{
  DTEvent tmp(at(0).signal, at(0).tag);
  std::vector<int> id_sig(m_sigType.size()), 
    ids_sig(m_sigType.size()), 
    id_tag(m_tagType.size()),
    ids_tag(m_tagType.size());

  TTree* outputTree = new TTree(name.c_str(),name.c_str());
  for(size_t i = 0 ; i < m_sigType.size(); ++i )
  {
    outputTree->Branch((particleName(m_sigType, i)+"_PX").c_str(), &tmp.signal[4*i+0]); 
    outputTree->Branch((particleName(m_sigType, i)+"_PY").c_str(), &tmp.signal[4*i+1]); 
    outputTree->Branch((particleName(m_sigType, i)+"_PZ").c_str(), &tmp.signal[4*i+2]); 
    outputTree->Branch((particleName(m_sigType, i)+"_E").c_str(),  &tmp.signal[4*i+3]);
    outputTree->Branch((particleName(m_sigType, i)+"_ID").c_str(), &id_sig[i]);
    ids_sig[i] = ParticlePropertiesList::get( m_sigType[i] )->pdgID();
  }  
  for(size_t i = 0 ; i < m_tagType.size(); ++i )
  {
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PX").c_str(), &tmp.tag[4*i+0]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PY").c_str(), &tmp.tag[4*i+1]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_PZ").c_str(), &tmp.tag[4*i+2]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_E").c_str(),  &tmp.tag[4*i+3]); 
    outputTree->Branch(("Tag_"+particleName(m_tagType, i)+"_ID").c_str(), &id_tag[i]);
    ids_tag[i] = ParticlePropertiesList::get( m_tagType[i] )->pdgID();
  }
  for( auto& evt: *this ){
    bool swap = gRandom->Uniform() > 0.5;
    tmp.set(evt.signal, evt.tag);
    if( swap ) tmp.invertParity();
    for(size_t i=0; i != m_sigType.size(); ++i)
      id_sig[i] = swap ? -ids_sig[i] : ids_sig[i];
    for(size_t i=0; i != m_tagType.size(); ++i)
      id_tag[i] = swap ? -ids_tag[i] : ids_tag[i];
    outputTree->Fill();
  }
  return outputTree;
}    