#ifndef AMPGEN_LINESHAPES_H
#define AMPGEN_LINESHAPES_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Tensor.h"
/**
  @defgroup Lineshapes Lineshapes
  @brief Lineshapes are semi-empirical complex functions for describing the propagation and decay of short-lived resonances. 
  The lineshapes are stored in a static factory and can be accessed by key from anywhere in the program.
  For example, consider a @f$K^{*}(892)^{0}@f$ resonance decaying into a kaon and a pion, the propagation of which could be described 
  using the Relativistic Breit Wigner (BW) lineshape. The corresponding decay descriptor would be

  \code{cpp}
  K*(892)[BW]{K+,pi-}
  \endcode

  where in this case, the BW is redundant, as this the default if no alternative is specified. Now consider the broad scalar @f$K^{*}(1430)^{0}@f$ 
  meson, which cannot readily be described by a simple Breit Wigner due to a quasi nonresonant amplitude (alternatively, due to the presence of an 
  additional very broad scalar at lower masses, the @f$\kappa(800)@f$). 
  One alternative parameterisation of this system is the LASS lineshape 
  which can be used as

  \code{cpp}
  K(0)*(1430)[LASS]{K+,pi-}
  \endcode

  Lineshapes for describing different systems are detailed below. 
  Some, such as the LASS shape or kMatrix, pertain to the specific particles and final states. 
  Others, such as CoupledChannel or Spline can flexibly describe a variety of systems, but therefore require a more detailed specification by the user. 
  
  Generally, the relativistic Breit-Wigner a sufficient description for many purposes for nonscalar resonances that decay to a single two-body final state or are relatively far from relevant thresholds.More complex systems require more specialised descriptions. 
  */

#define DECLARE_LINESHAPE( X )                                                                                  \
  class X : public AmpGen::ILineshape {                                                                         \
    static std::string _id;                                                                                     \
    public:                                                                                                     \
    X(){ DEBUG("Constructing lineshape") ;}                                                                     \
    AmpGen::Expression get(const AmpGen::Expression& s, const AmpGen::Expression& s1,                           \
                                   const AmpGen::Expression& s2, const std::string& particleName,               \
                                   const unsigned int& L, const std::string& lineshapeModifier,                 \
                                   AmpGen::DebugSymbols* dbexpressions = 0) const override;                     \
    AmpGen::Expression get(const AmpGen::Expression& s, const std::vector<AmpGen::Tensor>& p,                   \
                           const std::string& particleName,                                                     \
                           const unsigned int& L, const std::string& lineshapeModifier,                         \
                           AmpGen::DebugSymbols* dbexpressions = nullptr ) const override;                      \
  }

#define DEFINE_LINESHAPE( X )                                                                                   \
  REGISTER_WITH_KEY( ILineshape, Lineshape::X, #X, std::string );                                               \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const std::vector<AmpGen::Tensor>& p,      \
                                        const std::string& particleName,                                        \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const { return                    \
                                           get(s, dot(p[0],p[0]), dot(p[1],p[1]),                               \
                                           particleName, L, lineshapeModifier, dbexpressions) ;}                 \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const AmpGen::Expression& s1,              \
                                        const AmpGen::Expression& s2, const std::string& particleName,          \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const

#define DEFINE_GENERIC_SHAPE( X )                                                                               \
  REGISTER_WITH_KEY( ILineshape, Lineshape::X, #X, std::string );                                               \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const AmpGen::Expression& s1,              \
                                        const AmpGen::Expression& s2, const std::string& particleName,          \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const { return 0;}                \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const std::vector<AmpGen::Tensor>& p,      \
                                        const std::string& particleName,                                        \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const 

namespace AmpGen
{
  class  ILineshape {
    public:
    virtual ~ILineshape() = default;
    virtual Expression get( const Expression& s, const Expression& s1, const Expression& s2,
                            const std::string& particleName, const unsigned int& L,
                            const std::string& lineshapeModifier, DebugSymbols* dbexpressions = nullptr ) const = 0;

    virtual Expression get( const Expression& s, const std::vector<AmpGen::Tensor>& p, const std::string& particleName,
                            const unsigned int& L, const std::string& lineshapeModifier,
                            AmpGen::DebugSymbols* dbexpressions = nullptr ) const = 0; 
    ILineshape* create() { return this; }
  };

  /** @ingroup Lineshapes namespace Lineshape 
      Namespace that contains all lineshapes, i.e. propagators for describing amplitudes and phases for resonances (and nonresonant) contributions to a total amplitude. 
   */

   namespace Lineshape
  {
    class LineshapeFactory : public AmpGen::Factory<ILineshape>
    {
    public:
      static Expression getLineshape( const std::string& lineshape, const Expression& s, const Expression& s1,
                                      const Expression& s2, const std::string& particleName, const unsigned int& L,
                                      std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
      static Expression getGenericShape(const std::string& lineshape,
                                        const Expression& s,  
                                        const std::vector<Tensor>& p, const std::string& particleName,
                                        const unsigned int& L, AmpGen::DebugSymbols* dbexpressions );
      static bool isLineshape( const std::string& lineshape );
    };
    
    DECLARE_LINESHAPE( None );
    
    /** @ingroup Lineshapes class BW 
        @brief Simple two-body Breit-Wigner lineshape that describes relatively narrow, isolated resonances that couple to a single channel / orbital configuration.
        
        The propagator as a function of the invariant-mass squared @f$s@f$ for a two-body final state with relative orbital angular momentum @f$l@f$ is given by
        @f[
            \mathcal{A}(s) = \frac{ k(m,\Gamma_0) B_L(q(s),0) }{ m^2 - s - i m_0 \Gamma_{l}(s) },
        @f] <br> 
        <B> Parameters: </B> <br>  
        <DD>
             @f$m@f$        : <EM>particleName_</EM>mass   Breit-Wigner mass, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)  <br>
             @f$\Gamma_0@f$ : <EM>particleName_</EM>width   Breit-Wigner width, defined as the width of resonance at the Breit-Wigner mass <br>
             @f$r@f$        : <EM>particleName_</EM>radius   Hadronic radius for Blatt-Weisskopf form-factor (defaults to 1.5GeV for light resonances, 3.5GeV for charm) <br>
        </DD> <br>
        <B> Modifiers: </B> <br>
        <DD>
        <EM> BL </EM> : Use Blatt-Weisskopf factors normalised at @f$ \sqrt{s}=m @f$ (by default, normalised at @f$\sqrt{s}=0@f$)
        </DD>
        \image html figs/BW_combined.png "Modulus and phase of the Relativistic Breit-Wigner propagator, for @f$l={0,4}@f$, using the mass and nominal width of the @f$\rho@f$ meson" 
    */
    DECLARE_LINESHAPE( BW );
 
    /// Breit-Wigner lineshape with fixed width
    DECLARE_LINESHAPE( SBW );
    /// Non-relativistic Breit-Wigner lineshape
    DECLARE_LINESHAPE( NonRelBW );

    /** @ingroup Lineshapes class GounarisSakurai
        @brief Gounaris-Sakurai lineshape models dispersive corrections to the I=1,J=1 \f$\pi\pi\f$ propagator (G.J. Gounaris and J.J. Sakurai, Phys. Rev. Lett.21, 244 (1968))
        
        <B> Parameters: </B> <br>  
        <DD>
             @f$m@f$        : <EM>particleName_</EM>mass   Breit-Wigner mass, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)  <br>
             @f$\Gamma_0@f$ : <EM>particleName_</EM>width   Breit-Wigner width, defined as the width of resonance at the Breit-Wigner mass <br>
             @f$r@f$        : <EM>particleName_</EM>radius   Hadronic radius for Blatt-Weisskopf form-factor (defaults to 1.5GeV for light resonances, 3.5GeV for charm) <br>
        </DD> <br>
        \image html figs/GS_combined.png "Gounaris-Sakurai lineshape for the ϱ(770) meson, with the equivalent relativistic Breit Wigner lineshape shown for comparison.
    */
    DECLARE_LINESHAPE( GounarisSakurai );

    /// Description of the \f$ K\pi \f$ S-wave, based on the fits to scattering data.
    DECLARE_LINESHAPE( LASS );

    /// Lineshape to describe resonances with coupled channels such as \f$f_{0}(980)^{0} / a_{0}(980) \f$ (S.M.Flatté, Phys. Lett B. 63, 224 (1976))
    DECLARE_LINESHAPE( Flatte );
    DECLARE_LINESHAPE( Bugg );
    DECLARE_LINESHAPE( Isotensor );

    /// Exponential form-factor of type \f$ e^{ - q^2 R^2 } \f$
    DECLARE_LINESHAPE( ExpFF );

    DECLARE_LINESHAPE( FormFactor );

    /// Gaussian lineshape \f$ = e^{ -(x-\mu)^2 / 2\sigma^{2} } \f$
    DECLARE_LINESHAPE( Gaussian );

    /** @ingroup Lineshapes class Poly
     *  @brief Polynominal shape \f$ \mathcal{A}(s) = \sum^n_i c_i s^{i} \f$ where the sum is to lineshapeModifier::Degree, and the free parameters of the shape are lineshapeModifier_ci 
     */
    DECLARE_LINESHAPE( Poly );

    /** @ingroup Lineshapes class kMatrix 
        @brief Anisovich-Sarantsev Isoscalar K-matrix from https://arxiv.org/abs/hep-ph/0204328

        Describes the isoscalar @f$ \pi\pi, KK, 4\pi \eta\eta, \eta\eta^\prime@f$ S-wave in terms of a five-by-five K-matrix and corresponding P-vector couplings.
        Includes a large number of parameters that can be fixed from the above publication. 
        These parameters can be found in the options directory, which in turn can be includes in the fit by adding 

        \code{cpp}
          Import $AMPGENROOT/options/kMatrix.opt
        \endcode 
        
        to the user configuration file. 
     */ 
    DECLARE_LINESHAPE( kMatrix );

    /** @ingroup Lineshapes class FOCUS
     *  @brief K matrix amplitudes used for I=1/2 and I=3/2 in the description of the \f$ K\pi \f$ S-wave in the analysis of @f$ D^{+}\rightarrow K^{-}\pi^{+}\pi^{+}@f$ https://arxiv.org/abs/0705.2248
     */
     DECLARE_LINESHAPE( FOCUS );

    /// I=1/2 and I=3/2 K Matrices used in the description of the \f$ K\pi\f$ S-wave in the analysis of @f$\eta_{c}\rightarrow K^{+} K^{-} \pi^{0}@f$ https://arxiv.org/abs/1701.04881
    DECLARE_LINESHAPE( PALANO );

    /** @ingroup Lineshapes class ObelixRho
     *  @brief Amplitude to describe the vector-isovector system, otherwise known as the @f$ \rho @f$ mesons. WARNING untested. 

         Vector-Isovector amplitude @f$(I=1, J=1)@f$ using a K-matrix to describe the @f$\pi\pi,  KK, \pi\pi\pi\pi @f$ channels using three poles, commonly associated with 
     the @f$ \rho(770), \rho(1450), \rho(1900) @f$ resonances. 
     */
    DECLARE_LINESHAPE( ObelixRho );

    /// K matrix to describe \f$K_1(1270) / K_1(1400)\f$. WARNING incompleted. 
    DECLARE_LINESHAPE( AxialKaon );
    
    /** @ingroup Lineshapes class kMatrixSimple 
        @brief Simple and flexible K matrix that implements a variable number of scalar channels and poles.   
        Flexible K matrix implementation that accepts a variable number of scalar channels and pole terms.
        Generally only likely to be useful for pedagoical examples. */
    DECLARE_LINESHAPE( kMatrixSimple );
    
    /** @ingroup Lineshapes class MIPWA
        @brief Model-Independent Partial Wave parameterisation using cubic splines.
      
        The real and imaginary (or amplitude and phase) of the lineshape are fixed at \f$ N \f$ discrete positions to independent pairs of free parameters to be determined in a fit, 
        and the value of the lineshape determined elsewhere by interpolating between these values using cubic splines. 
        Based on techniques first used by the E791 collaboration in studying @f$ K \pi @f$ scalars, see https://arxiv.org/abs/hep-ex/0510045. <br>
        The amplitude is given by 
         @f[
            \mathcal{A}(s) = \sum_n ( a_n + b_n \tilde{s} + c_n \tilde{s}^2 + d_n \tilde{s}^3 ) \mathrm{rect}\left( \frac{s - nL}{nL} \right)
         @f]
        <B> Parameters: </B> <br>
      <DD>
         @f$\mathcal{R}_j@f$ : <EM>particleName</EM>::Spline::Re::j  The real value of the lineshape evaluated at the @f$ j @f$th knot. <br>
         @f$\mathcal{I}_j@f$ : <EM>particleName</EM>::Spline::Im::j  The imaginary value of the lineshape evaluated at the @f$ j @f$th knot.  <br>
      </DD> <br>    
     */
    DECLARE_LINESHAPE( MIPWA );
    /** @ingroup Lineshapes class GSpline 
        @brief 
      */
    
    DECLARE_LINESHAPE( GSpline );
    DECLARE_LINESHAPE( FormFactorSpline );
    DECLARE_LINESHAPE( DecaySpline );
    DECLARE_LINESHAPE( InelasticSpline );

    /** @ingroup Lineshapes class CoupledChannel 
        @brief Description of a resonance that decays to multiple two and three-body final states. 
      */
    DECLARE_LINESHAPE( CoupledChannel );
    /** @ingroup Lineshapes class GenericKmatrix
        @brief Implementation of a generic K-matrix
      */
    DECLARE_LINESHAPE(GenericKmatrix);
    /// ``Lineshape'' that implements the Dalitz plot distribution for decays \f$ \eta \rightarrow \pi^{+}\pi^{-}\pi^{0}\f$
    DECLARE_LINESHAPE( EtaDalitz );

  } // namespace Lineshape
  
  Expression Q2( const Expression& Msq, const Expression& M1sq, const Expression& M2sq );

  Expression kFactor( const Expression& mass, const Expression& width,
                      std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
  Expression BlattWeisskopf_Norm( const Expression& z, const Expression& z0, unsigned int L );

  Expression BL( const Expression& s, const Expression& s0, const Expression& s1, const Expression& s2,
                 const Expression& radius, const unsigned int& L );

  Expression BlattWeisskopf( const Expression& z, unsigned int L );

  Expression pol( const AmpGen::Expression& X, const std::vector<Expression>& p );

  std::vector<Expression> parameterVector( const std::string& name, const unsigned int& nParam );

  Expression width( const Expression& s, const Expression& s1, const Expression& s2, const Expression& mass,
                    const Expression& width, const Expression& radius, unsigned int L,
                    DebugSymbols* dbexpressions = nullptr );

} // namespace AmpGen
#endif
