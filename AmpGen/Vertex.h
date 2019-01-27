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

/** \defgroup Vertices Spin
  * Vertices for the decays of particles are formed from the polarisation vectors or spinors of the decay products, 
  * and the momentum and momentum transfer of the decaying particle. 
  * These objects are combined using a vertex to form a generalised current that can be used to compute the angular momentum components of the transition matrix.  
  * In general, these objects are only described for decays of the type \f$ 1\to 2 \f$, with more complex structures constructed from such quasi two-body decays (see Isobar model), 
  * but more general amplitudes can also be defined. 
  * \author T.Evans
  */

/** \ingroup Vertices \macro DECLARE_VERTEX
  * Macro to declare a vertex 
  */
#define DECLARE_VERTEX( NAME )                                                                                         \
  struct NAME : public VertexBase {                                                                                       \
    NAME(){ DEBUG("Constructing vertex");                    }                                                         \
    virtual Tensor operator()( const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2, DebugSymbols* db = 0 ) override;  \
    static std::string _id;                                                                                            \
  }

#define DEFINE_VERTEX( VERTEX )                                                                                   \
  REGISTER_WITH_KEY( Vertex::VertexBase, Vertex::VERTEX, #VERTEX, std::string );                                                     \
  Tensor Vertex::VERTEX::operator()( const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2, DebugSymbols* db )

namespace AmpGen
{
  namespace Vertex
  {
    /** \ingroup Vertices class VertexBase
     *  Virtual base class from which all the other vertices derive, in essence this is 
     *  just a named function pointer that can create a pointer to itself, i.e. such 
     *  that it can be constructed using \ref Factory "Factory".
     */
    struct VertexBase {
      
      /* Calculate the generalised current for this decay process, as a function of:
      * \param P The momentum of the decaying particle
      * \param Q The momentum transfer between the two decay products
      * \param V1 The polarisation tensor or spinor of the first decay product. By construction, the 
      * particle ordering convention used implies that the first decay product has the higher spin
      * \param V2 The polarisation tensor or spinor of the second decay product. 
      * \param db Optional collection of debug symbols for evaluating this amplitude. 
      */
      virtual AmpGen::Tensor operator()(const AmpGen::Tensor& P, 
                                        const AmpGen::Tensor& Q, 
                                        const AmpGen::Tensor& V1, 
                                        const AmpGen::Tensor& V2,
                                        AmpGen::DebugSymbols* db = nullptr ) = 0;

      virtual ~VertexBase() = default;
      VertexBase* create() { return this; }   
    };
    /// \ingroup Vertices class S_SS_S 
    /// \brief \f$ S_1 S_2 \f$
    DECLARE_VERTEX( S_SS_S );

    /// \ingroup Vertices class S_VV_S 
    /// \brief \f$ g_{\mu\nu} V_1^\mu V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_S );
    
    /// \ingroup Vertices class S_VV_S1
    /// \brief \f$ S_{\mu\nu} V_1^\mu V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_S1 );

    /// \ingroup Vertices class S_VV_P
    /// \f$ \varepsilon_{\alpha\beta\mu\nu} P^{\alpha} L^{\beta} V_1^{\mu} V_2^{\nu} \f$
    DECLARE_VERTEX( S_VV_P );
    
    /// \ingroup Vertices class S_VV_D
    /// \brief \f$ L_{\mu\nu} V_1^\mu V_2^\nu \f$
    DECLARE_VERTEX( S_VV_D );

    /// \ingroup Vertices class S_VS_P
    /// \brief \f$ L_{\mu} V_1^{\mu} S_2 \f$
    DECLARE_VERTEX( S_VS_P );

    /// \ingroup Vertices class S_TV_P
    /// \brief \f$ L^{\mu} T_{\mu\nu} V^{\nu} \f$
    DECLARE_VERTEX( S_TV_P );

    /// \ingroup Vertices class S_TV_D 
    /// \brief \f$ \varepsilon_{\mu\nu\alpha\beta} T^{\mu\gamma} L_{\gamma}^{\nu} P^{\alpha} V^{\beta} \f$
    DECLARE_VERTEX( S_TV_D );
    
    /// \ingroup Vertices class S_TS_D 
    /// \brief \f$ T^{\mu\nu} L_{\mu\nu}\f$
    DECLARE_VERTEX( S_TS_D );
    DECLARE_VERTEX( S_TT_S );

    DECLARE_VERTEX( V_SS_P );
    DECLARE_VERTEX( V_VS_P );
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
    DECLARE_VERTEX( V_ff_P );
    DECLARE_VERTEX( V_ff_P1 );  
  } // namespace Vertex
  class VertexFactory : public Factory<Vertex::VertexBase>
  {
  public:
    static Tensor getSpinFactor( const Tensor& P, const Tensor& Q, const Tensor& V1, const Tensor& V2,
                                 const std::string& name, DebugSymbols* db = nullptr );
    static Tensor getSpinFactorNBody( const std::vector<std::pair<Tensor, Tensor>>& tensors, const unsigned int& mL,
                                      DebugSymbols* db = nullptr );
    static bool isVertex( const std::string& hash );
  };

  /// \ingroup Vertices function Orbital_PWave 
  /// Helper function that computes the L=1 orbital momentum operator, i.e. 
  /// \f$ L_{\mu} = q_{\mu} - \frac{p_{\nu}q^{\nu} p_{\mu} }{ p_{\alpha} p^{\alpha}} \f
  Tensor Orbital_PWave( const Tensor& A, const Tensor& B );
  Tensor Orbital_DWave( const Tensor& A, const Tensor& B );
  Tensor Spin1Projector( const Tensor& A );
  Tensor Spin2Projector( const Tensor& A );

  Tensor Spin1hProjector( const Tensor& B );
  Tensor Spin3hProjector( const Tensor& A );

  Tensor gamma_twiddle( const Tensor& P );
  Tensor Gamma4Vec();
  Tensor Bar( const Tensor& P );
  Tensor slash( const Tensor& P );
} // namespace AmpGen
#endif
