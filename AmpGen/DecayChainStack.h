#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stack>

#include "AmpGen/EventType.h"
#include "AmpGen/Particle.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/SmallVector.h"
#include "AmpGen/Event.h"
#include <TRandom3.h>

namespace AmpGen { 

  class DecayChainStackBase
  {
    public:
      DecayChainStackBase() = default;
      DecayChainStackBase( const Particle& particle ); 
      virtual ~DecayChainStackBase(){};
      virtual Event makeEvent(TRandom3* rndm) const = 0; 
      virtual double maxWeight() const = 0; 
      virtual unsigned NP()      const = 0; 
      virtual bool operator==( const DecayChainStackBase& other ) const = 0;  
  };

  template <unsigned N> class DecayChainStack : public DecayChainStackBase {

    struct Node {
      enum Type { BW, Flat, Stable, QuasiStable }; 
      Type type       = Type::Flat; 
      unsigned l;
      unsigned r; 
      std::pair<double, double> range; 
      std::pair<double, double> phi_range; 
      double bwMass   = 0; 
      double bwWidth  = 0;
      SmallVector<unsigned, N> lfs; 
      SmallVector<unsigned, N> rfs; 
      Node(Type type = Type::Flat, const unsigned& l=0, const unsigned& r=0) : type(type), l(l), r(r) {};
      void print() const 
      {
        std::string type_string = "";
        switch( type )
        {
          case Type::BW          : type_string = "BW         "; break; 
          case Type::Flat        : type_string = "Flat       "; break; 
          case Type::Stable      : type_string = "Stable     "; break; 
          case Type::QuasiStable : type_string = "QuasiStable"; break; 
        };
        INFO( type_string << " decay indices: [" << l << ", " << r << "], range : [" << range.first << ", " << range.second << "] "  
            ", [" << vectorToString( lfs, ",") << "], [" << vectorToString( rfs, ",") << "]" ); 
      }
      template <typename T>
        T operator()( const T& s ) const 
        {
          switch( type ) {
            case Type::BW          : 
              return  ( bwMass* bwWidth ) /( (phi_range.second-phi_range.first )*  ( (s - bwMass*bwMass)*(s-bwMass*bwMass) + bwMass*bwMass*bwWidth*bwWidth) );

            case Type::Flat        : return 1/(range.second - range.first); 
            case Type::Stable      : return 1; 
            case Type::QuasiStable : return 1; 
          };
          return 0; 
        }

      double s_nom() const 
      {
        switch( type ) {
          case Type::BW          : return bwMass * bwMass; 
          case Type::Flat        : return (range.second + range.first ) * 0.5; 
          case Type::Stable      : return range.first; 
          case Type::QuasiStable : return range.first; 
        };
        return 0; 
      }
      double operator()( TRandom3* random ) const 
      {
        switch( type ) {
          case Type::BW :
            {
              auto y = tan( ( phi_range.second - phi_range.first ) * random->Rndm() ); 
              auto mG = bwMass * bwWidth; 
              auto s0 = bwMass * bwMass; 
              return s0 + mG * ( mG * y + range.first - s0 ) / ( mG - y * ( range.first - s0 ) );  
            }
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
        INFO("Setting particle: " << particle->name() << " " << particle->isStable() << " " << particle->isQuasiStable() << " " << particle->lineshapeContains({"FormFactor", "None","TD"})  );  
        if( particle->isStable() ) type = Type::Stable; 
        else if( particle->isQuasiStable() ) type = Type::QuasiStable; 
        else if( particle->lineshapeContains({"FormFactor", "None","TD"}) ) type == Type::Flat; 
        else { 
          type = Type::BW; 
          phi_range.first  = atan((range.first  - bwMass*bwMass)/(bwMass*bwWidth));
          phi_range.second = atan((range.second - bwMass*bwMass)/(bwMass*bwWidth));
        }
      }
      void setRange( const double& min, const double& max ) 
      {
        if( type == Type::BW or type == Type::Flat ) 
        {
          range = std::make_pair(min, max);
          if( type == Type::BW )
          {
            phi_range.first  = atan((range.first  - bwMass*bwMass)/(bwMass*bwWidth));
            phi_range.second = atan((range.second - bwMass*bwMass)/(bwMass*bwWidth));
          }
        }
      }
      bool operator==( const Node& other ) const 
      {
        return  type  == other.type && 
          range == other.range &&
          phi_range == other.phi_range &&  
          bwMass == other.bwMass && 
          bwWidth == other.bwWidth &&
          l == other.l && 
          r == other.r && 
          lfs == other.lfs && 
          rfs == other.rfs; 
      }
    }; /// each node sets the energy / momentum of a set of particles 
    std::array<double, N> m_m0; 
    std::array<Node, N-1> m_nodes;
    double m_rhoMax=1;

    public:

    DecayChainStack( const Particle& particle )
    {
      auto fs = particle.getFinalStateParticles(); 
      double minMass = particle.mass(); 
      for( int i = 0 ; i != N; ++i ){ m_m0[i] = fs[i]->mass(); minMass -= m_m0[i]; }
      std::stack<std::pair<std::shared_ptr<Particle>, unsigned>> toDo; 
      toDo.emplace( std::make_shared<Particle>(particle), 0); 
      unsigned counter = 0; 
      while( toDo.size() != 0 )
      {
        auto [current,index] = toDo.top();
        toDo.pop(); 
        if( current->daughters().size() != 2 )
        {  
          std::vector<Particle> particles;
          auto dp = current->daughters(); 
          for( unsigned int i = 1; i != dp.size(); ++i ) particles.push_back( *dp[i] );
          current = std::make_shared<Particle>( current->name(), *dp[0], Particle("NonResS0",particles) );   
        }
        INFO( index << " " << *current ); 
        if( current->daughters().size() != 2 ) FATAL("Should be fixed..");
        auto fs  = current->quasiStableTree().daughters(); 
        auto min_mass = std::accumulate(fs.begin(), fs.end(), 0., [](double acc, auto& p){ return acc + p->mass(); } );
        auto G = current->isQuasiStable() ? 0 : current->props()->mass() * current->props()->width() * 10; 
        auto s0 = current->mass() * current->mass();
        auto d1_index = current->daughter(0)->index();
        auto d2_index = current->daughter(1)->index();
        if( !current->daughter(0)->isStable() ){ d1_index = ++counter + N; toDo.emplace( current->daughter(0), d1_index-N); } 
        if( !current->daughter(1)->isStable() ){ d2_index = ++counter + N; toDo.emplace( current->daughter(1), d2_index-N); }   
        m_nodes[index] = Node(Node::Type::Flat, d1_index, d2_index);
        if( index == 0 or current->isQuasiStable() ) m_nodes[index].range = std::make_pair(s0 - G, s0 + G ); 
        else m_nodes[index].range = std::make_pair( pow(min_mass,2), pow(minMass + min_mass,2 ) );  
        m_nodes[index].set( &(*current) ); 
      }
      if( NamedParameter<bool>("DecayChainStack::AggressiveOptimisation", false )  )
      {
        WARNING("Using aggressive optimisation of phase space, can cause problems for relative normalisation of different decay topologies"); 
        /// this optimisation 'cuts' the phase space of the a decay at the maximal extent of its sibling, i.e. reduces 
        /// the physically allowed region somewhat. 
        /// To be checked: this probably changes the normalisation of different topologies so shouldn't be switched on by default 
        /// But probably works for single topology cases ... 
        for( auto& n : m_nodes ) 
        {
          if( n.l > N && n.r > N ){
            m_nodes[n.l-N].setRange( m_nodes[n.l-N].range.first, n.range.second + m_nodes[n.r-N].range.first - 2 * std::sqrt( n.range.second * m_nodes[n.r-N].range.first ) ); 
            m_nodes[n.r-N].setRange( m_nodes[n.r-N].range.first, n.range.second + m_nodes[n.l-N].range.first - 2 * std::sqrt( n.range.second * m_nodes[n.l-N].range.first ) ); 
          }
        } 
      }
      auto concat = [](auto v1, const auto& v2 ){ v1.insert( v1.end(), v2.begin(), v2.end() ); return v1; }; 
      for(auto n = m_nodes.rbegin(); n != m_nodes.rend(); ++n ){
        auto s1 = n->l < N ? m_m0[n->l]*m_m0[n->l] : m_nodes[n->l - N].range.first; 
        auto s2 = n->r < N ? m_m0[n->r]*m_m0[n->r] : m_nodes[n->r - N].range.first;

        auto rho = rho_2(n->range.second, s1, s2 );
        // INFO(m_nodes.size() - ( n-m_nodes.rbegin()) + N -1  << " " << rho << " " << n->l << " " << n->r << " " << n->range.first << " " << n->range.second << " " << n->type <<
        //     " s(max): " << n->range.second << 
        //     " s1(min): " << s1  << " s2(min): " << s2 );
        m_rhoMax *= rho; 
        n->lfs   = n->l < N ? SmallVector<unsigned, N>{n->l} : concat( m_nodes[n->l - N].lfs, m_nodes[n->l - N].rfs);
        n->rfs   = n->r < N ? SmallVector<unsigned, N>{n->r} : concat( m_nodes[n->r - N].lfs, m_nodes[n->r - N].rfs);
      }
      for( auto& node : m_nodes ) node.print();
      if( m_rhoMax <= 0 ) ERROR("RhoMax^2 < 0! " << m_rhoMax );
      m_rhoMax = std::sqrt(m_rhoMax);
    }
    virtual double maxWeight() const { return m_rhoMax; }
    virtual unsigned NP() const { return N; } 
    virtual bool operator==( const DecayChainStackBase& other ) const 
    {
      auto op = dynamic_cast<const DecayChainStack<N>*>( &other );
      return op && this->m_m0 == op->m_m0 && this->m_nodes == op->m_nodes; 
    }
    /*
       double genPdf(const double* event) const
       {
       double g = 1; 
       for( const auto& node : m_nodes ){
    //INFO( "s(" << vectorToString( node.lfs, ",") << "," << vectorToString(node.rfs,",") << ") = " << s(event, node.lfs, node.rfs) ); 
    g *= node(s(event, node.lfs, node.rfs));
    }
    return g;    
    }
    */
    template <typename T>
      T genPdf(const T* event) const
      {
        T g = 1; 
        for( const auto& node : m_nodes ){
          //INFO( "s(" << vectorToString( node.lfs, ",") << "," << vectorToString(node.rfs,",") << ") = " << s(event, node.lfs, node.rfs) ); 
          g *= node(s(event, node.lfs, node.rfs));
        }
        return g;    
      }
    template <typename T> 
      T s(const T* event, 
          const SmallVector<unsigned, N>& lfs, 
          const SmallVector<unsigned, N>& rfs) const 
      {
        std::array<T,4> P = {0.}; 
        for( const auto& it : lfs )
        {
          P[0] += event[4*it + 0 ];
          P[1] += event[4*it + 1 ];
          P[2] += event[4*it + 2 ];
          P[3] += event[4*it + 3 ];
        }
        for( const auto& it : rfs )
        {
          P[0] += event[4*it + 0 ];
          P[1] += event[4*it + 1 ];
          P[2] += event[4*it + 2 ];
          P[3] += event[4*it + 3 ];
        }
        return P[3] * P[3] - P[0] * P[0] - P[1] * P[1] - P[2] * P[2];  
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
    void set( double* output, double px, double py, double pz, double E) const 
    {
      output[0] = px; 
      output[1] = py; 
      output[2] = pz; 
      output[3] = E; 
    }

    std::pair<double, std::array<double, 2*N-1>> proposal( TRandom3* rndm ) const 
    {
      std::pair<double, std::array<double, 2*N-1>> rt; 
      auto& [weight, state ] = rt; 
      for( unsigned int i = 0 ; i != N ; ++i ) state[i] = m_m0[i] * m_m0[i]; 
      weight = 1; 
      for( int i = m_nodes.size()-1 ; i >= 0; --i ){
        state[i + N]   = m_nodes[i](rndm);
        const auto& s  = state[i+N];
        const auto& sl = state[m_nodes[i].l];
        const auto& sr = state[m_nodes[i].r];
        double v = rho_2(s, sl, sr) * ( s > sl + sr + 2 * std::sqrt(sl*sr) ); 
        weight  = v < 0 ? 0 : v * weight; 
      }
      weight = std::sqrt( weight ); 
      return rt; 
    }
    void debug( TRandom3* rndm ) const 
    {
      std::pair<double, std::array<double, 2*N-1>> rt; 
      auto& [weight, state ] = rt; 
      for( unsigned int i = 0 ; i != N ; ++i ) state[i] = m_m0[i] * m_m0[i]; 
      weight = 1; 
      for( int i = m_nodes.size()-1 ; i >= 0; --i ){
        state[i + N]   = m_nodes[i](rndm);
        const auto& s  = state[i+N];
        const auto& sl = state[m_nodes[i].l];
        const auto& sr = state[m_nodes[i].r];
        double v = rho_2(s, sl, sr) * ( s > sl + sr + 2 * std::sqrt(sl*sr) ); 
        INFO( "Node[" << i << "] w = " << v << "s[" << m_nodes[i].l << " " << m_nodes[i].r << "] = " << s << " [ sl=  " << sl << " sr= " << sr << "]" ); 
        weight  = v < 0 ? 0 : v * weight; 
      }
      weight = std::sqrt( weight ); 
    }

    void fill_from_state( double* event, const std::array<double, 2*N-1>& state, TRandom3* rndm ) const 
    {
      for( int i = 0 ; i != N; ++i ) event[4*i+3] = m_m0[i];
      for(auto n = m_nodes.rbegin(); n != m_nodes.rend(); ++n)
      {
        const auto& sl = state[n->l];
        const auto& sr = state[n->r]; 
        const auto& s  = state[ &(*n) - &(*m_nodes.begin()) + N ];
        const auto nz  = 2 *  rndm->Rndm() - 1;
        const auto phi = 2*M_PI*rndm->Rndm();      
        const auto p   = 0.5 * std::sqrt( s * rho_2(s, sl, sr ) ); 
        const auto sZ  = std::sqrt( 1  - nz*nz ); 
        const auto nx  = sZ * std::cos(phi); 
        const auto ny  = sZ * std::sin(phi); 
        const auto gl  = std::sqrt( 1 + p*p/sl);  
        const auto vgl =  p / std::sqrt(sl); 
        const auto gr  = std::sqrt( 1 + p*p/sr);  
        const auto vgr =  p / std::sqrt(sr); 

        if( n->lfs.size == 1 ) set( event + 4 * n->lfs[0], p*nx, p*ny, p*nz, std::sqrt( sl + p*p ) );  
        else for( const auto& l : n->lfs ) boost(event + 4 *l,  nx , ny,  nz, gl, vgl );

        if( n->rfs.size == 1 ) set( event + 4 *n->rfs[0], -p*nx, -p*ny, -p*nz, std::sqrt( sr + p*p ) );        
        else for( const auto& l : n->rfs ) boost(event + 4 *l, -nx, -ny, -nz, gr, vgr ); 
      }
    }

    virtual AmpGen::Event makeEvent(TRandom3* rndm) const
    {
      std::array<double, 2 * N - 1> state = {0}; 
      double rho_v = 1 ; 
      do {
        auto rt = proposal( rndm ); 
        rho_v = std::get<0>(rt);
        state = std::get<1>(rt); 
      } while( rho_v < m_rhoMax * rndm->Uniform() ); 
      AmpGen::Event event( N * 4 ); 
      fill_from_state( event, state, rndm ); 
      return event; 
    }
  };

  DecayChainStackBase* make_decay_chain_stack( const Particle& particle )
  {
    auto fs = particle.getFinalStateParticles();
    switch ( fs.size() )
    { 
      case(1 ) : return new DecayChainStack<1>( particle ); 
      case(2 ) : return new DecayChainStack<2>( particle ); 
      case(3 ) : return new DecayChainStack<3>( particle ); 
      case(4 ) : return new DecayChainStack<4>( particle ); 
      case(5 ) : return new DecayChainStack<5>( particle ); 
      case(6 ) : return new DecayChainStack<6>( particle ); 
      case(7 ) : return new DecayChainStack<7>( particle ); 
      case(8 ) : return new DecayChainStack<8>( particle ); 
      case(9 ) : return new DecayChainStack<9>( particle ); 
      case(10) : return new DecayChainStack<10>( particle ); 
    }
    return nullptr; 
  }

}
