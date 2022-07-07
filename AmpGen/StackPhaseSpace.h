#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stack>

#include "AmpGen/EventType.h"
#include "AmpGen/Particle.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/EventList.h"

namespace AmpGen { 

  class DecayChainStackBase
  {
    public:
      DecayChainStackBase() = default;
      DecayChainStackBase( const Particle& particle ); 
      virtual Event makeEvent(TRandom3* rndm) const = 0; 
  };


  template <typename type, unsigned max_size> 
    struct small_vector {
      typedef type value_type; 

      std::array<type, max_size> store; 
      unsigned size = 0; 
      auto begin() const { return store.begin(); }
      auto end()   const { return store.begin() + size; } 
      auto begin()  { return store.begin(); }
      auto end()    { return store.begin() + size; } 
      value_type& operator[]       (unsigned i)       { return store[i]; } 
      const value_type& operator[] (unsigned i) const { return store[i]; } 
      void push_back( const value_type& thing ){ store[size++] = thing; }
      template <typename iterator_type, typename other_iterator_type> 
        void insert(iterator_type mbegin, other_iterator_type ibegin, other_iterator_type iend )
        {  
          for( auto it = ibegin; it != iend; ++it ){
            *( mbegin + (it - ibegin) ) = *it;  
          }    
          size += (iend - ibegin ); 
        }
      small_vector() = default; 
      small_vector( std::initializer_list<value_type>&& values )
      {
        insert( begin(), values.begin(), values.end() );
      }
    };


  template <unsigned N> class DecayChainStack : public DecayChainStackBase {

    struct Node {
      enum Type { BW, Flat, Stable, QuasiStable }; 
      unsigned l;
      unsigned r; 
      std::pair<double, double> range; 
      double bwMass   = 0; 
      double bwWidth  = 0; 
      double phiMin   = 0; 
      double phiMax   = 0; 
      Type type       = Type::Flat; 
      Node(Type type = Type::Flat, const unsigned& l=0, const unsigned& r=0) : type(type), l(l), r(r) {};
      small_vector<unsigned, N> lfs; 
      small_vector<unsigned, N> rfs; 
      double s_nom() const 
      {
        switch( type ) {
          case Type::BW : return bwMass * bwMass; 
          case Type::Flat : return (range.second + range.first ) * 0.5; 
          case Type::Stable : return range.first; 
          case Type::QuasiStable : return range.first; 
        };
        return 0; 
      }
      double operator()( TRandom3* random ) const 
      {
        switch( type ) {
          case Type::BW :
            return bwMass * bwMass + bwMass * bwWidth * tan( (phiMax - phiMin ) * random->Rndm() +  phiMin );
          case Type::Flat : 
            return (range.second-range.first) * random->Rndm() + range.first;
          case Type::Stable:
            return range.first;   
          case Type::QuasiStable:
            return range.first;   
        } 
        return 0; 
      }
      void set( const Particle* particle ) 
      {
        auto props = particle->props();
        bwWidth = props->width();
        bwMass = props->mass();
        if( particle->isStable() ) type = Type::Stable; 
        else if( particle->isQuasiStable() ) type = Type::QuasiStable; 
        else if( particle->lineshape().find("BW") != std::string::npos ){
          type = Type::BW; 
          phiMin = atan((range.first  - bwMass*bwMass)/(bwMass*bwWidth));
          phiMax = atan((range.second - bwMass*bwMass)/(bwMass*bwWidth));
        }
        else type = Type::Flat;
      }
    }; /// each node sets the energy / momentum of a set of particles 
    std::array<double,N>  m_m0; 
    std::array<Node, N-1> m_nodes;
    double m_rhoMax=1;

    public: 

    DecayChainStack( const Particle& particle )
    {
      auto fs = particle.getFinalStateParticles(); 
      double minMass = particle.mass(); 
      for( int i = 0 ; i != N; ++i ){ m_m0[i] = fs[i]->mass(); minMass -= m_m0[i]; }
      std::stack<const Particle*> toDo; 
      toDo.push(&particle); 
      unsigned int counter = N; 
      while( toDo.size() != 0 )
      {
        auto current = toDo.top();
        toDo.pop(); 
        auto& d1 = *current->daughter(0);
        auto& d2 = *current->daughter(1); 
        auto fs = current->getFinalStateParticles();
        auto min_mass = std::accumulate(fs.begin(), fs.end(), 0., [](double acc, auto& p){ return acc + p->mass(); } );
        auto G = current->isQuasiStable() ? 0 : current->props()->mass() * current->props()->width() * 10;
        m_nodes[counter - N]     = Node(Node::Type::Flat, d1.isStable() ? d1.index() : counter+1, d2.isStable() ? d2.index() : counter+1+ !d1.isStable());
        if( counter == N ) m_nodes[counter -N].range = std::make_pair( current->mass()*current->mass() - G, current->mass() * current->mass() + G ); 
        else if( current->isQuasiStable() ) m_nodes[counter -N].range = std::make_pair( current->mass() * current->mass(), current->mass() * current->mass() );
        else m_nodes[counter-N].range = std::make_pair( pow(min_mass,2), pow(minMass + min_mass,2 ) );  
        if( current->name() == "D*(2010)+" or current->name() == "D*(2007)0" ) 
          m_nodes[counter -N].range = std::make_pair( current->mass() * current->mass() - G, current->mass() * current->mass() + G );
        m_nodes[counter-N].set( current ); 
        if( !d1.isStable() ) toDo.push( &d1 ); 
        if( !d2.isStable() ) toDo.push( &d2 ); 
        counter++; 
      }
      auto concat = [](auto v1, const auto& v2 ){ v1.insert( v1.end(), v2.begin(), v2.end() ); return v1; }; 

      for(auto n = m_nodes.rbegin(); n != m_nodes.rend(); ++n ){
        auto rho = rho_2(n->range.second, n->l < N ? m_m0[n->l]*m_m0[n->l] : m_nodes[n->l].range.first, n->r < N ? m_m0[n->r]*m_m0[n->r] : m_nodes[n->r].range.first );       
        m_rhoMax *= rho; 
        n->lfs       = n->l < N ? small_vector<unsigned, N>{n->l} : concat( m_nodes[n->l - N].lfs, m_nodes[n->l - N].rfs);
        n->rfs       = n->r < N ? small_vector<unsigned, N>{n->r} : concat( m_nodes[n->r - N].lfs, m_nodes[n->r - N].rfs);
      }
      m_rhoMax = sqrt(m_rhoMax);
    }

    double rho_2( const double& s0, const double& s1, const double& s2 ) const 
    {
      return 1 - 2 * ( s1 + s2 )/s0 + ( s1 - s2 )*( s1 - s2 ) / (s0*s0); 
    }
    void boost( double* output, double nx, double ny, double nz, double g, double vg ) const
    {
      double nv = output[0] * nx + output[1] * ny + output[2] * nz;
      output[0] += ( (g-1) * nv + vg * output[3] ) * nx ;
      output[1] += ( (g-1) * nv + vg * output[3] ) * ny ;
      output[2] += ( (g-1) * nv + vg * output[3] ) * nz ;
      output[3] = g * output[3]  + vg * nv;
    }
    virtual AmpGen::Event makeEvent(TRandom3* rndm) const
    {
      std::array<double, 2 * N - 1> state; 
      for( unsigned int i = 0 ; i != N ; ++i ) state[i] = m_m0[i] * m_m0[i]; 
      double rho_v = 1 ; 
      do {
        rho_v = 1; 
        for( int i = m_nodes.size()-1 ; i >= 0; --i ){
          state[i + N] = m_nodes[i](rndm); 
          double v = rho_2(state[i+N], state[m_nodes[i].l], state[m_nodes[i].r] ); 
          rho_v  = v < 0 ? 0 : v * rho_v; 
        }
      } while( sqrt(rho_v) < m_rhoMax * rndm->Uniform() ); 
      AmpGen::Event event( N * 4 ); 
      for( int i = 0 ; i != N; ++i ) event[4*i+3] = m_m0[i];
      for(auto n = m_nodes.rbegin(); n != m_nodes.rend(); ++n)
      {
        const auto& sl = state[n->l];
        const auto& sr = state[n->r]; 
        const auto& s  = state[ &(*n) - &(*m_nodes.begin()) + N ];
        const auto nz  = 2 *  rndm->Rndm() - 1;
        const auto phi = 2*M_PI*rndm->Rndm();      
        const auto p   = 0.5 * sqrt( s * rho_2(s, sl, sr ) ); 
        const auto sZ  = std::sqrt( 1  - nz*nz ); 
        const auto nx  = sZ * std::cos(phi); 
        const auto ny  = sZ * std::sin(phi); 
        const auto gl  = sqrt( 1 + p*p/sl);  
        const auto vgl =  p / sqrt(sl); 
        const auto gr  = sqrt( 1 + p*p/sr);  
        const auto vgr =  p / sqrt(sr); 
        for( const auto& l : n->lfs ) boost(event + 4 *l,  nx , ny,  nz, gl, vgl );
        for( const auto& l : n->rfs ) boost(event + 4 *l, -nx, -ny, -nz, gr, vgr ); 
      }
      return event; 
    }
  };
}
