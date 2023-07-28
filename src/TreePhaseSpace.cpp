#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Event.h"
#include "AmpGen/DecayChainStack.h"
#include "AmpGen/simd/utils.h"
#include "TRandom3.h"
#include <numeric>
#include <random>

using namespace AmpGen;

TreePhaseSpace::TreePhaseSpace(const Particle& decayChain, const EventType& type, TRandom* rndm) :
  m_rand(rndm == nullptr ? (TRandom3*)gRandom : (TRandom3*)rndm),
  m_type(type)
{
  auto orderings = decayChain.identicalDaughterOrderings();
  for( auto& ordering : orderings )
  {
    Particle p = decayChain;
    p.setOrdering(ordering);
    m_gen.push_back(make_decay_chain_stack(p));
  }
  initialise_weights(); 
  setRandom(m_rand);
}

TreePhaseSpace::TreePhaseSpace(const std::vector<Particle>& decayChains, const EventType& type, TRandom* rndm) :
  m_rand(rndm == nullptr ? (TRandom3*)gRandom : (TRandom3*)rndm),
  m_type(type)
{
  auto add_decay_chain = [this]( Particle p ) mutable 
  {
    auto orderings = p.identicalDaughterOrderings();
    DEBUG("Adding tree: " << p );
    for( auto& ordering : orderings )
    {
      p.setOrdering(ordering);
      auto s = make_decay_chain_stack(p);
      bool any_of_fast = std::any_of( this->m_gen.begin(), this->m_gen.end(), [s](const auto it){ return *it == *s; } ); 
      if(any_of_fast) {continue;  }
      this->m_gen.push_back(s); 
    }
  }; 
  for( const auto& decayChain : decayChains ) add_decay_chain( decayChain );  
//   add_decay_chain( Particle( decayChains[0].quasiStableTree() ) );
  
  initialise_weights(); 
  setRandom(rndm);
  DEBUG("w(max) = " << m_wmax << "#chains: " << m_gen.size() );  
}

void TreePhaseSpace::initialise_weights() 
{
  auto weights = std::vector<double>( m_gen.size(), 1./m_gen.size() );
  m_wmax = (*std::max_element( m_gen.begin(), m_gen.end(), [](auto a, auto b){ return a->maxWeight() > b->maxWeight();} )) -> maxWeight(); 
  m_dice = DiscreteDistribution(weights); 
}

template <unsigned N> void fill_event(Event& event, 
                                      const std::vector<DecayChainStackBase*>& phsp, 
                                      TRandom3* rndm, 
                                      const DiscreteDistribution& weight_distribution, 
                                      double wmax ) 
{
  std::pair<double, std::array<double, 2 * N -1>> rt;
  const DecayChainStack<N>* node = nullptr;
  if( dynamic_cast<const DecayChainStack<N>*>( phsp[0] ) == 0 )
    FATAL("Failed to cast: " << phsp[0]->NP() << " " << N );

  unsigned nTrial = 0; 
  do {
    node = static_cast<const DecayChainStack<N>*>(phsp[weight_distribution(rndm)]); 
    rt = node->proposal( rndm ); 
    nTrial++;
    if( nTrial > 1000 ) {
      ERROR("Tried " << nTrial << " and still no dice (pun intentional)" << " " << std::get<0>(rt) ); 
      node->debug(rndm); 
    }
  } while( std::get<0>(rt) < rndm->Uniform() * wmax ); 
  node->fill_from_state(event.address(), std::get<1>(rt), rndm );

  double genPdf = 0; 
  for( unsigned int i = 0 ; i != phsp.size(); ++i )
  {
    genPdf += weight_distribution.probabilities()[i] * 
              static_cast<const DecayChainStack<N>*>(phsp[i])->genPdf( event.address() );
  }  
  event.setGenPdf( genPdf );
//  event.print(); 
}

template <unsigned N>
real_v get_gen_pdf( const real_v* value, const std::vector<DecayChainStackBase*>& phsp, const DiscreteDistribution& weight_distribution )
{
  real_v genPdf = 0;  
  for( unsigned j = 0 ; j != phsp.size(); ++j )
    genPdf += weight_distribution.probabilities()[j] * static_cast<const DecayChainStack<N>*>(phsp[j])->genPdf(value); 
  return genPdf; 
}

template <unsigned N> std::vector<double> calculate_weights(const EventListSIMD& events, 
                                                            const std::vector<DecayChainStackBase*>& phsp, 
                                                            const DiscreteDistribution& weight_distribution)
{
  // procedure proposed in: https://arxiv.org/pdf/hep-ph/9405257.pdf
  std::vector<real_v> weights( phsp.size(), 0 );
  for( unsigned block = 0 ; block < events.nBlocks(); ++block )
  {
    auto genPdf = get_gen_pdf<N>( events.block(block), phsp, weight_distribution );    
    for( unsigned j = 0 ; j != phsp.size(); ++j )
      weights[j] += static_cast<const DecayChainStack<N>*>( phsp[j] )->genPdf( events.block(block)) * events.genPDF(block) * events.genPDF(block) / genPdf; 
  }
  std::vector<double> return_weights( phsp.size() );
  for(unsigned i = 0 ; i != phsp.size(); ++i )
    return_weights[i] = utils::sum_elements( weights[i] );
  return return_weights; 
}

Event TreePhaseSpace::makeEvent()
{
  auto s= m_type.size(); 
  Event event( s * 4 );
  switch (s)
  {
    case(1) :  fill_event<1>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(2) :  fill_event<2>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(3) :  fill_event<3>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(4) :  fill_event<4>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(5) :  fill_event<5>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(6) :  fill_event<6>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(7) :  fill_event<7>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(8) :  fill_event<8>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(9) :  fill_event<9>( event, m_gen, m_rand, m_dice, m_wmax ); break; 
    case(10) : fill_event<10>(event, m_gen, m_rand, m_dice, m_wmax ); break; 
  }
  return event; 
}

void TreePhaseSpace::recalculate_weights(const EventListSIMD& events)
{
  std::vector<double> weights; 
  switch( m_type.size() )
  {
    case(1)  : weights = calculate_weights<1>( events, m_gen, m_dice); break; 
    case(2)  : weights = calculate_weights<2>( events, m_gen, m_dice); break;
    case(3)  : weights = calculate_weights<3>( events, m_gen, m_dice); break;
    case(4)  : weights = calculate_weights<4>( events, m_gen, m_dice); break;
    case(5)  : weights = calculate_weights<5>( events, m_gen, m_dice); break;
    case(6)  : weights = calculate_weights<6>( events, m_gen, m_dice); break;
    case(7)  : weights = calculate_weights<7>( events, m_gen, m_dice); break;
    case(8)  : weights = calculate_weights<8>( events, m_gen, m_dice); break; 
    case(9)  : weights = calculate_weights<9>( events, m_gen, m_dice);  break;
    case(10) : weights = calculate_weights<10>( events, m_gen, m_dice); break;
  }
  for( int j = 0 ; j != m_gen.size(); ++j )
  {
    weights[j] = pow( weights[j] / events.size(), 0.5)  * m_dice.probabilities()[j]; 
  }
  m_dice = DiscreteDistribution(weights); 
}

void TreePhaseSpace::setRandom( TRandom* rand )
{
  m_rand = (TRandom3*)rand;
}

EventType TreePhaseSpace::eventType() const 
{
  return m_type; 
}

size_t TreePhaseSpace::size() const 
{
  return m_gen.size();
}

TreePhaseSpace::~TreePhaseSpace() 
{
  for( auto generator : m_gen ) delete generator; 
}
