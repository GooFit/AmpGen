#ifndef AMPGEN_IVERTEX_H
#define AMPGEN_IVERTEX_H

#include <map>
#include <utility>
#include <vector>
#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/MsgService.h"

/** @defgroup Vertices Spin
    Vertices for the decays of particles are formed from the polarisation vectors or spinors of the decay products, 
    and the momentum and momentum transfer of the decaying particle. 
    These objects are combined using a vertex to form a generalised current that can be used to compute the angular momentum components of the transition matrix.  
    In general, these objects are only described for decays of the type \f$ 1\to 2 \f$, with more complex structures constructed from such quasi two-body decays (see Isobar model), 
    but more general amplitudes can also be defined. 
  */

/** @ingroup Vertices macro DECLARE_VERTEX
  * Macro to declare a vertex 
  */
#define DECLARE_VERTEX(NAME)                                                                                         \
  struct NAME : public Base {                                                                                       \
    NAME(){ DEBUG("Constructing vertex");                    }                                                         \
    virtual Tensor operator()(const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2, DebugSymbols* db = 0 ) override;  \
    static std::string _id;                                                                                            \
  }

#define DEFINE_VERTEX(VERTEX)                                                                                   \
  REGISTER_WITH_KEY( Vertex::Base, Vertex::VERTEX, #VERTEX, std::string );                                                     \
  Tensor Vertex::VERTEX::operator()( const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2, DebugSymbols* db )

namespace AmpGen
{
  /** @ingroup Vertices namespace Vertex 
      Namespace that contains the base class for vertices, Vertex::Base, as well as the implementations 
      of specific spin couplings and some helper functions such as the orbital operators. 
    */
  namespace Vertex
  {
    /** @ingroup Vertices class Base
        @brief Base class for all spin vertices.  
      Virtual base class from which all the other vertices derive, in essence this is 
       just a named function pointer that can create a pointer to itself, i.e. such 
       that it can be constructed using a Factory.
     */
    struct Base {
      
      /** Calculate the generalised current for this decay process, as a function of:
         @param P The momentum of the decaying particle
         @param Q The momentum transfer between the two decay products
         @param V1 The polarisation tensor or spinor of the first decay product. By construction, the 
         particle ordering convention used implies that the first decay product has the higher spin
         @param V2 The polarisation tensor or spinor of the second decay product. 
         @param db Optional collection of debug symbols for evaluating this amplitude. 
      */
      virtual AmpGen::Tensor operator()(const AmpGen::Tensor& P, 
                                        const AmpGen::Tensor& Q, 
                                        const AmpGen::Tensor& V1, 
                                        const AmpGen::Tensor& V2,
                                        AmpGen::DebugSymbols* db = nullptr ) = 0;

      virtual ~Base() = default;
      Base* create() { return this; }   
    };
    /// \ingroup Vertices class S_SS_S 
    /// \brief \f$ S = S_1 S_2 \f$
    DECLARE_VERTEX( S_SS_S );

    /// \ingroup Vertices class S_VV_S 
    /// \brief \f$ S = g_{\mu\nu} V_1^\mu V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_S );
    
    /// \ingroup Vertices class S_VV_S1
    /// \brief \f$ S = S_{\mu\nu} V_1^\mu V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_S1 );

    /// \ingroup Vertices class S_VV_P
    /// \f$ S = \varepsilon_{\alpha\beta\mu\nu} P^{\alpha} L^{\beta} V_1^{\mu} V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_P );
    
    /// \ingroup Vertices class S_VV_D
    /// \brief \f$ S = L_{\mu\nu} V_1^\mu V_2^\nu \f$
    DECLARE_VERTEX( S_VV_D );

    /// \ingroup Vertices class S_VS_P
    /// \brief \f$ S = L_{\mu} V_1^{\mu} S_2 \f$
    DECLARE_VERTEX( S_VS_P );

    /// \ingroup Vertices class S_TV_P
    /// \brief \f$ S = L^{\mu} T_{\mu\nu} V^{\nu} \f$
    DECLARE_VERTEX( S_TV_P );

    /// \ingroup Vertices class S_TV_D 
    /// \brief \f$ S = \varepsilon_{\mu\nu\alpha\beta} T^{\mu\gamma} L_{\gamma}^{\nu} P^{\alpha} V^{\beta} \f$
    DECLARE_VERTEX( S_TV_D );
    
    /// \ingroup Vertices class S_TS_D 
    /// \brief \f$ S = T^{\mu\nu} L_{\mu\nu}\f$
    DECLARE_VERTEX( S_TS_D );
    
    /// \ingroup Vertices class S_TT_S 
    /// \brief \f$ S = T_1^{\mu\nu} T_{2\mu\nu}\f$
    DECLARE_VERTEX( S_TT_S );

    /// @ingroup Vertices class V_SS_P
    /// @brief @f$ V^{\mu} = L^{\mu} S_1 S_2 @f$
    DECLARE_VERTEX( V_SS_P );

    DECLARE_VERTEX( V_VS_P );
    
    /// @ingroup Vertices class V_SS_P
    /// @brief @f$ V^{\mu} = S^{\mu\nu} V_{1\nu} S_2 @f$
    DECLARE_VERTEX( V_VS_S );
    DECLARE_VERTEX( V_VS_D );
    DECLARE_VERTEX( V_TS_P );
    DECLARE_VERTEX( V_TS_D );

    DECLARE_VERTEX( T_VS_D );
    DECLARE_VERTEX( T_VS_P );
    DECLARE_VERTEX( T_SS_D );
    DECLARE_VERTEX( T_TS_D );
    DECLARE_VERTEX( T_TS_S );

    DECLARE_VERTEX( f_fS_S );
    DECLARE_VERTEX( f_fS_S1 );
    DECLARE_VERTEX( f_fS_P );
    DECLARE_VERTEX( f_fS_P1 );

    DECLARE_VERTEX( f_Vf_S );
    DECLARE_VERTEX( f_Vf_S1 );
    DECLARE_VERTEX( f_Vf_P );
    DECLARE_VERTEX( f_Vf_P1 );
    DECLARE_VERTEX( f_Vf_P2 );
    DECLARE_VERTEX( f_Vf_P3 );
    
    DECLARE_VERTEX( f_Vf_D );
    DECLARE_VERTEX( f_Vf_D1 );
    
    DECLARE_VERTEX( f_Tf_P );
    
    DECLARE_VERTEX( r_fS_P );
    DECLARE_VERTEX( r_fS_D );
    DECLARE_VERTEX( f_rS_P );
    DECLARE_VERTEX( f_rS_D );
    DECLARE_VERTEX( f_rS_P1 );

    DECLARE_VERTEX( S_ff_S );
    DECLARE_VERTEX( S_ff_S1 );
    DECLARE_VERTEX( V_ff_S );
    DECLARE_VERTEX( V_ff_S1 );  
    class Factory : public AmpGen::Factory<Vertex::Base>
    {
    public:
      static Tensor getSpinFactor( const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2,
                                   const std::string& name, DebugSymbols* db = nullptr );
      static Tensor getSpinFactorNBody( const std::vector<std::pair<Tensor, Tensor>>& tensors, const unsigned int& mL,
                                        DebugSymbols* db = nullptr );
      static bool isVertex( const std::string& hash );
    };
  } // namespace Vertex

  /** @ingroup Vertices function Orbital_PWave 
      @brief Helper function that computes the @f$L=1@f$ orbital momentum operator.
      Helper function that computes the @f$L=1@f$ orbital momentum operator, which is given by 
      @f[ L_{\mu} = q_{\mu} - p_{\mu} \frac{p_{\nu}q^{\nu}}{ p^2 }, @f]
      where @f$ p @f$ is total momentum and @f$ q @f$ is the momentum difference between 
      the two particles in the state.
   */
  Tensor Orbital_PWave(const Tensor& p, const Tensor& q);
  
  /** @ingroup Vertices function Orbital_DWave 
      @brief Helper function that computes the @f$L=2@f$ orbital momentum operator.  
      Helper function that computes the @f$L=2@f$ orbital momentum operator, which is given by
      @f[ L_{\mu\nu} = L_{\mu}L_{\nu} - \frac{L^2}{3} S_{\mu\nu}, @f]
      where @f$ L_{\mu} @f$ is the Orbital_PWave operator and @f$ S_{\mu\nu} @f$ is the 
      Spin1Projector. */
  Tensor Orbital_DWave(const Tensor& p, const Tensor& q);

  /** @ingroup Vertices function Spin1Projector
      @brief Helper function that computes the projection operator onto a spin one state.
      Helper function that projects some lorentz object onto the momentum @f$p_\mu @f$ of state, 
      and is given by 
      @f[ S_{\mu\nu} = g_{\mu\nu} - \frac{p_{\mu}p_{\nu}}{p^2}. @f] */
  Tensor Spin1Projector(const Tensor& p);
  
  /** @ingroup Vertices function Spin2Projector
      @brief Helper function that computes the projection operator onto a spin one state.
      Helper function that projects some lorentz object onto the momentum @f$p_\mu @f$ of state, 
      and is given by 
      @f[ S_{\mu\nu\alpha\beta} = \frac{1}{2}\left( S_{\mu\alpha} S_{\nu\beta} + S_{\mu\beta}S_{\nu\alpha} \right) - \frac{1}{3} S_{\mu\nu} S_{\alpha\beta}, @f]
      where @f$ S_{\mu\nu} @f$ is the spin-one projection operator (see Spin1Projector). */
  Tensor Spin2Projector(const Tensor& p);
  
  /** @ingroup Vertices function Spin1hProjector
      @brief Helper function that projects a spinor.
      Helper function that projects out a spin-half state 
      @f[ S_{ab} = \frac{1}{2m}\left( {p\!\!\!/} + m I \right) @f]
    */
  Tensor Spin1hProjector( const Tensor& B );
  Tensor Spin3hProjector( const Tensor& A );

  Tensor gamma_twiddle( const Tensor& P );
  Tensor Gamma4Vec();
  Tensor Bar( const Tensor& P );
  Tensor slash( const Tensor& P );
} // namespace AmpGen
#endif
