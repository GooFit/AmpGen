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
#include "AmpGen/Particle.h"

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
  class X : public AmpGen::Lineshape::Base {                                                                    \
    static std::string _id;                                                                                     \
    public:                                                                                                     \
    X(){ DEBUG("Constructing lineshape") ;}                                                                     \
    AmpGen::Expression get(const AmpGen::Expression& s, const AmpGen::Expression& s1,                           \
                                   const AmpGen::Expression& s2, const std::string& particleName,               \
                                   const unsigned int& L, const std::string& lineshapeModifier,                 \
                                   AmpGen::DebugSymbols* dbexpressions = 0) const override;                     \
    AmpGen::Expression get(const AmpGen::Particle& p,                                                           \
                           const std::string& lineshapeModifier,                                                \
                           AmpGen::DebugSymbols* dbexpressions = nullptr ) const override;                      \
  }

#define DEFINE_LINESHAPE( X )                                                                                   \
  REGISTER_WITH_KEY( Lineshape::Base, Lineshape::X, #X, std::string );                                          \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Particle& p,                                              \
                                        const std::string& lineshapeModifier,                                   \
                                        AmpGen::DebugSymbols* dbexpressions ) const { return                    \
                                           get(p.massSq(), p.daughter(0)->massSq(), p.daughter(1)->massSq(),    \
                                           p.name(), p.L(), lineshapeModifier, dbexpressions) ;}                \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const AmpGen::Expression& s1,              \
                                        const AmpGen::Expression& s2, const std::string& particleName,          \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const

#define DEFINE_GENERIC_SHAPE( X )                                                                               \
  REGISTER_WITH_KEY( Lineshape::Base, Lineshape::X, #X, std::string );                                          \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const AmpGen::Expression& s1,              \
                                        const AmpGen::Expression& s2, const std::string& particleName,          \
                                        const unsigned int& L, const std::string& lineshapeModifier,            \
                                        AmpGen::DebugSymbols* dbexpressions ) const { return 0;}                \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Particle& p,                                              \
                                        const std::string& lineshapeModifier,                                   \
                                        AmpGen::DebugSymbols* dbexpressions ) const 

namespace AmpGen
{
  class Particle; // forward definitions

  /** @ingroup Lineshapes namespace Lineshape 
      Namespace that contains all lineshapes, i.e. propagators for describing amplitudes and phases for resonances (and nonresonant) contributions to a total amplitude. 
   */

  namespace Lineshape
  {
    class  Base {
      public:
      virtual ~Base() = default;
      virtual Expression get( const Expression& s, const Expression& s1, const Expression& s2,
                              const std::string& particleName, const unsigned& L,
                              const std::string& lineshapeModifier, DebugSymbols* dbexpressions = nullptr ) const = 0;
      virtual Expression get( const AmpGen::Particle& p, 
                              const std::string& lineshapeModifier,
                              AmpGen::DebugSymbols* dbexpressions = nullptr ) const = 0; 
      Base* create() { return this; }
    };

    class Factory : public AmpGen::Factory<Lineshape::Base>
    {
    public:
      static Expression get(const std::string& lineshape, const Expression& s, const Expression& s1,
                            const Expression& s2, const std::string& particleName, const unsigned& L,
                            std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
      static Expression get(const std::string& lineshape,
                            const Particle& p,
                            AmpGen::DebugSymbols* dbexpressions );
      static bool isLineshape( const std::string& lineshape );
    };
    
    DECLARE_LINESHAPE( None );
    
    /** @ingroup Lineshapes class BW 
        @brief Simple two-body Breit-Wigner lineshape that describes relatively narrow, isolated resonances that couple to a single channel / orbital configuration.
        
        The propagator as a function of the invariant-mass squared @f$s@f$ for a two-body final state with relative orbital angular momentum @f$l@f$ is given by
        @f[
            \mathcal{A}(s) = \frac{ k(m,\Gamma_0) B_l(qr,0) }{ m^2 - s - i m \Gamma_{l}(s) },
        @f] 
        where the width @f$ \Gamma_{l}(s) @f$ is given by
        @f[
            \Gamma_{l}(s) = \frac{\Gamma_0 q m_0}{q_0 \sqrt{s}} \left(\frac{q}{q_0}\right)^{2l} B_l( qr, q_0 r )^2,
        @f] 
        where @f$q@f$ is the linear momentum of either decay product in the rest frame of the parent, while @f$q_0@f$ the same quantity evaluated at the resonance mass. 
        The Blatt-Weisskopf functions, @f$B_l(q,q_0)@f$ control the behaviour of functions at large momentum transfers. The width is normalised such that @f$\Gamma(m^2) = \Gamma_0 @f$.  
        The parameters of the lineshape are tabulated below
         Parameter              | User name                            | Description 
         -----------------------|--------------------------------------|------------------------------------------------------------------------
         @f$m@f$                | <EM>particleName_</EM>mass           | Breit-Wigner mass, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)
         @f$\Gamma_0@f$         | <EM>particleName_</EM>width          | Breit-Wigner width, defined as the width of resonance at the Breit-Wigner mass
         @f$r@f$                | <EM>particleName_</EM>radius         | Hadronic radius for Blatt-Weisskopf form-factor (defaults to 1.5GeV for light resonances, 3.5GeV for charm)
        
        <EM> BL </EM> : Use Blatt-Weisskopf factors normalised at @f$ \sqrt{s}=m @f$ (by default, normalised at @f$\sqrt{s}=0@f$)
        \image html BW_combined.png "Modulus and phase of the Relativistic Breit-Wigner propagator, for orbital momentum up-to four, using the mass and nominal width of the rho meson." 
    */
    DECLARE_LINESHAPE( BW );
 
    /** @ingroup Lineshapes class SBW 
        @brief Breit-Wigner lineshape with fixed width
        
        Basic Breit-Wigner lineshape with a fixed width and without a form factor.
        @f[
            \mathcal{A}(s) = \frac{ k(m,\Gamma_0) }{ m^2 - s - i m \Gamma_0 },
        @f] 
        The parameters of the lineshape are tabulated below
         Parameter              | User name                            | Description 
         -----------------------|--------------------------------------|------------------------------------------------------------------------
         @f$m@f$                | <EM>particleName_</EM>mass           | Breit-Wigner mass, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)
         @f$\Gamma_0@f$         | <EM>particleName_</EM>width          | Breit-Wigner width, defined as the width of resonance at the Breit-Wigner mass
     */
    DECLARE_LINESHAPE( SBW );
    
    /** @ingroup Lineshapes class NonRelBW 
        @brief Nonrelativistic Breit-Wigner lineshape 

       Non-relativistic Breit-Wigner lineshape, given by 

      @f[
       \mathcal{A}(s) = \sqrt{ \frac{ \Gamma_0 }{ 2 \pi } } \left( \sqrt{s} - m_0 - i \Gamma_0 /2 \right)^{-1}
      @f]
    */  
    DECLARE_LINESHAPE( NonRelBW );

    /** @ingroup Lineshapes class GounarisSakurai
        @brief Gounaris-Sakurai lineshape models dispersive corrections to the @f$I=1,J=1@f$ @f$\pi\pi@f$ propagator.
        
        The Gounaris-Sakurai lineshape models dispersive corrections to the I=1,J=1 @f$\pi\pi@f$ propagator (G.J. Gounaris and J.J. Sakurai, Phys. Rev. Lett.21, 244 (1968))
        @f[
            \mathcal{A}(s) = \frac{ k(m,\Gamma_0) B_l(qr,0) }{ m^2(s) - s - i m \Gamma(s) },
        @f]
        where the running mass @f$ m^2(s)@f$ is given by
        @f[
          m^2(s) = m^2 + \frac{\Gamma_0}{\pi} \left( \frac{ 2 q^3 m }{q_0^3 \sqrt{s} } \log\left(\frac{\sqrt{s} + q }{ 2 m_\pi }\right) +
           \frac{ q^2 ( s - 3m^2) + s(m^2-s) }{ m q_0^2 } \log\left(\frac{m+q_0}{2 m_\pi}\right) +
           \frac{m^2 - s}{q_0}
          \right)
        @f]

        Parameter              | User name                            | Description 
        -----------------------|--------------------------------------|------------------------------------------------------------------------
        @f$m@f$                | <EM>particleName_</EM>mass           | Breit-Wigner mass, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)  <br>
        @f$\Gamma_0@f$         | <EM>particleName_</EM>width          | Breit-Wigner width, defined as the width of resonance at the Breit-Wigner mass <br>
        @f$r@f$                | <EM>particleName_</EM>radius         | Hadronic radius for Blatt-Weisskopf form-factor (defaults to 1.5GeV for light resonances, 3.5GeV for charm) <br>
        
         \image html figs/GS_combined.png "Gounaris-Sakurai lineshape for the ϱ(770) meson, with the equivalent relativistic Breit Wigner lineshape shown for comparison.
    */
    DECLARE_LINESHAPE( GounarisSakurai );

    /** @ingroup Lineshapes class LASS 
        @brief Description of the @f$ K\pi @f$ S-wave, based on the fits to scattering data.
        The LASS parameterisation of the @f$ K\pi @f$ S-wave is derived from fits to ~ elastic @f$ K \pi @f$ scattering data, which is approximately up to the @f$ K \eta^\prime @f$ threshold. 
        In this regime, unitarity implies that phases, rather than amplitudes should be summed. 
        In this context, a slow varying nonresonant phase,
         @f[
           \tan(\phi_{NR}) = \frac{2 a q}{2 + arq^2},
         @f]
         where @f$a, r@f$ are parameters determined from scattering data, is added to the phase shift of a Breit-Wigner, 
         @f[
           \tan(\phi_{BW}) = \frac{m\Gamma(s)}{m^2 -s },
         @f] 
         normally associated to the @f$ K^*(1430)^{0} @f$ resonance.
         The total amplitude is therefore:
         @f[
          \mathcal{A}(s)  = \frac{2 a \sqrt{s} }{ 2 + arq^2 - 2iaq} + \frac{2+arq^2 + 2iaq}{2+arq^2 - 2iaq }\mathcal{A}_{BW}(s) = \mathcal{A}_{NR}(s) + \mathcal{A}^{\prime}_{BW}(s).
         @f]
         As this expression somewhat resembles the sum of a Breit-Wigner with a slowly varying nonresonant component, the two parts are sometimes split apart with an additional production amplitude placed on one or the other. 
         These can be accessed separately using the modifiers LASS.BW and LASS.NR, and given independent coupling constants. From the user interface:
         \code{.cpp}
           K*(1430)bar0[LASS.NR]{K-,pi+}  0 2.0 0.1 0 -90  0.1
           K*(1430)bar0[LASS.BW]{K-,pi+}  0 1.5 0.1 0  20  0.1
         \endcode 
         Corresponds to the overall lineshape: 
         @f[
           \mathcal{A}(s)
           = \left(\texttt{K*(1430)bar0[LASS.NR]{K-,pi+}}\right) \mathcal{A}_{NR}(s) + 
             \left(\texttt{K*(1430)bar0[LASS.BW]{K-,pi+}}\right) \mathcal{A}^{\prime}_{BW}(s) 
           = 2.0 e^{-90^\mathrm{o}} \mathcal{A}_{NR}(s) + 1.5 e^{-20^\mathrm{o}} \mathcal{A}^{\prime}_{BW}(s).
         @f]
         An alternative way to introduce a different coupling between the two components is to fix one coupling and to introduce an additional production coupling. 
         In the user interface, this might for example in @f$ D^{0} \to K^{-} \pi^{+} \pi^{0} @f$ decays:
         \code{.cpp}
           D0{K*(1430)bar0,pi0}           0 1.3 0.0 0  15  0.1
           K*(1430)bar0[LASS.NR]{K-,pi+}  0 2.0 0.1 0 -90  0.1
           K*(1430)bar0[LASS.BW]{K-,pi+}  2 1.0 0.0 0   0  0.0
         \endcode 
         Then the lineshape with production couplings will be: 
         @f[
           \mathcal{A}(s)
           = \left(\texttt{D0{K*(1430)bar0,pi0}}\right)\left(
           = \left(\texttt{K*(1430)bar0[LASS.NR]{K-,pi+}}\right) \mathcal{A}_{NR}(s) + 
             \mathcal{A}^{\prime}_{BW}(s) \right)
           = 1.3 e^{15^\mathrm{o}} \left( 2.0 e^{-90^\mathrm{o}} \mathcal{A}_{NR}(s) + \mathcal{A}^{\prime}_{BW}(s) \right).
         @f]
        
        Parameter              | User name                            | Description 
        -----------------------|--------------------------------------|------------------------------------------------------------------------
        @f$m@f$                | <EM>particleName_</EM>mass           | Breit-Wigner mass of the resonant component, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)  <br>
        @f$\Gamma_0@f$         | <EM>particleName_</EM>width          | Breit-Wigner width of the resonant component, defined as the width of resonance at the Breit-Wigner mass <br>
        @f$a@f$                | LASS::a                              | Scattering length of the nonresonant component, defaults to @f$2.07\mathrm{G\kern -0.1em eV}^{-1}@f$
        @f$r@f$                | LASS::r                              | Effective interaction length of the nonresonant component, defaults to @f$3.32\mathrm{G\kern -0.1em eV}^{-1}@f$
    */  
    DECLARE_LINESHAPE( LASS );
    /** @ingroup Lineshapes class GLASS 
        @brief Generalised description of the @f$ K\pi@f$ S-wave, taken from BABAR / BELLE analysis of @f$ D^{0} \to K_{S}^{0} \pi^{+} \pi^{-} @f$ https://journals.aps.org/prd/pdf/10.1103/PhysRevD.98.112012 
        
        Alternative parameterisation of the @f$ K \pi @f$ S-wave taken from the BABAR/BELLE analysis of @f$ D^{0} \to K_{S}^{0} \pi^{+} \pi^{-} @f$ https://journals.aps.org/prd/pdf/10.1103/PhysRevD.98.112012. 
        It introduces four additional parameters, two amplitudes and two phases, and is written 
        @f[
          \mathcal{A}(s) = R\sin\delta_R e^{i(\delta_R + 2\delta_F)} + F\sin\delta_F e^{i\delta_F}
        @f]
        where 
        @f[
          \delta_R = \phi_R + \tan^{-1}\left(\frac{m\Gamma(s)}{m^2 -s }\right)
        @f]
        and 
        @f[
          \delta_F = \phi_F + \tan^{-1}\left(\frac{2 a q}{2 + arq^2}\right),
        @f]
        that is, the 'resonant' and 'nonresonant' phases of the LASS formulation with additional 'proudction' phases @f$ \phi_R @f$ and @f$\phi_F @f$, and production couplings @f$R, F@f$. One of the production couplings can typically be fixed to 1, normally the resonant coupling, @f$R @f$.  

        Parameter              | Interface parameter name             | Description 
        -----------------------|--------------------------------------|------------------------------------------------------------------------
        @f$m@f$                | <EM>particleName_</EM>mass           | Breit-Wigner mass of the resonant component, defined as energy at which the self-energy of the resonance is purely imaginary (defaults to value in PDG)  <br>
        @f$\Gamma_0@f$         | <EM>particleName_</EM>width          | Breit-Wigner width of the resonant component, defined as the width of resonance at the Breit-Wigner mass <br>
        @f$a@f$                | <EM>particleName</EM>::GLASS::a      | Scattering length of the nonresonant component, defaults to @f$2.07\mathrm{G\kern -0.1em eV}^{-1}@f$
        @f$r@f$                | <EM>particleName</EM>::GLASS::r      | Effective interaction length of the nonresonant component, defaults to @f$3.32\mathrm{G\kern -0.1em eV}^{-1}@f$
        @f$\phi_R@f$           | <EM>particleName</EM>::GLASS::phiR   | Additional production phase of the resonant component, defaults to @f$0^{\mathrm{o}}@f$, measured in degrees. 
        @f$\phi_F@f$           | <EM>particleName</EM>::GLASS::phiF   | Additional production phase of the nonresonant component, defaults to @f$0^{\mathrm{o}}@f$, measured in degrees. 
        @f$R@f$                | <EM>particleName</EM>::GLASS::R      | Additional production coupling of the resonant component, defaults to @f$1@f$.
        @f$F@f$                | <EM>particleName</EM>::GLASS::F      | Additional production coupling of the nonresonant component, defaults to @f$1@f$.
    */  
    DECLARE_LINESHAPE(GLASS);

    /** @ingroup Lineshapes class Flatte
        @brief Lineshape to describe resonances with coupled channels such as @f$f_{0}(980)^{0} / a_{0}(980) @f$.
        
        The lineshape was first described by S.M.Flatté in Phys. Lett B. 63, 224 (1976), and describes a single, isolated resonance that couples to a pair of 
        channels, which can have a large impact on the lineshape if the opening of one of the channels is near to the mass of the resonance. This lineshape 
        can be considered a special case of the K matrix or coupled channel formalisms, and is explicitly implemented with the couplings and channels for 
        the @f$f_{0}(980)^{0}@f$ and @f$a_{0}(980)^{0}@f$ resonances, which couple to the channels @f$\pi\pi@f$ and @f$KK@f$ for the @f$f_{0}(980)^{0}@f$ and @f$\pi\eta@f$ and @f$KK@f$ for the @f$a_{0}(980)^{0}@f$. In particular, the opening of the @f$KK@f$ threshold close to the resonance mass strongly distorts the lineshape. For a generic implementation of a coupled channel, see the CoupledChannel lineshape.  
        @f[
            \mathcal{A}(s) = \frac{1}{ m^2 - s - i m \Gamma(s) },
        @f]
        where the running width is given by 
       @f[
         \Gamma(s) = \frac{ g_{\pi\pi} }{s} \left( \Lambda^{1/2}(s,m_{\pi}^2,m_{\pi}^2)  + \frac{g_{KK}}{g_{\pi\pi}} \Lambda^{1/2}(s,m_K^2, m_K^2) \right)
       @f] 
       or 
       @f[
         \Gamma(s) = \frac{ g_{\pi\eta}}{s} \left( \Lambda^{1/2}(s,m_{\pi}^2,m_{\eta}^2)  + \frac{g_{KK}}{g_{\pi\eta}} \Lambda^{1/2}(s,m_K^2, m_K^2) \right) 
       @f]
       for the @f$f_0(980)^{0}@f$ and the @f$a_0(980)^{0}@f$, respectively. 

     */
    DECLARE_LINESHAPE( Flatte );
    DECLARE_LINESHAPE( Bugg );
    DECLARE_LINESHAPE( Isotensor );

    /// Exponential form-factor of type \f$ e^{ - q^2 R^2 } \f$
    DECLARE_LINESHAPE( ExpFF );

    DECLARE_LINESHAPE( FormFactor );

    /** @ingroup Lineshapes class Gaussian 
        @brief Gaussian shape for (relatively) long lived states that are limited by experimental resolution, rather than natural width.
        The gaussian lineshape has the form 
        @f[
          \mathcal{A}(s) = e^{ -(s-\mu)^2 / 2\sigma^2 },
        @f]
        which is potentially useful for fitting very long-lived contributions that are essentially limited by the experimental resolution, 
        rather than the natural lifetime. This type of parameterisation will assess contributions from interference incorrectly.

        Parameter              | Interface parameter name             | Description 
        -----------------------|--------------------------------------|------------------------------------------------------------------------
        @f$\mu@f$              | <EM>lineshapeModifier</EM>_mean      | Mean of the gaussian distribution.
        @f$\sigma@f$           | <EM>lineshapeModifier</EM>_sigma     | Width of the gaussian distribution. 
      */
    DECLARE_LINESHAPE( Gaussian );

    /** @ingroup Lineshapes class Poly
     *  @brief Polynominal shape \f$ \mathcal{A}(s) = \sum^n_i c_i s^{i} \f$ where the sum is to lineshapeModifier::Degree, and the free parameters of the shape are lineshapeModifier_ci 
     */
    DECLARE_LINESHAPE( PolyNR );
    
    /** @ingroup Lineshapes class PolyNR
     *  @brief Polynominal shape \f$ \mathcal{A}(s) = \sqrt{ \sum_{ij} c_{ij} s^{i} (s^\prime)^{j} } \f$ where the sum is to lineshapeModifier::Degree, and the free parameters of the shape are lineshapeModifier_cij 
     */
    DECLARE_LINESHAPE( Poly );

    /** @ingroup Lineshapes class FOCUS
     *  @brief K matrix amplitudes used for I=1/2 and I=3/2 in the description of the \f$ K\pi \f$ S-wave in the analysis of @f$ D^{+}\rightarrow K^{-}\pi^{+}\pi^{+}@f$ https://arxiv.org/abs/0705.2248
     */
     DECLARE_LINESHAPE( FOCUS );

    /// I=1/2 and I=3/2 K Matrices used in the description of the \f$ K\pi\f$ S-wave in the analysis of @f$\eta_{c}\rightarrow K^{+} K^{-} \pi^{0}@f$ https://arxiv.org/abs/1701.04881
    DECLARE_LINESHAPE( PALANO );

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
            \mathcal{A}(s) = \sum_n ( a_n + b_n \tilde{s} + c_n \tilde{s}^2 + d_n \tilde{s}^3 ) \mathrm{rect}\left( \frac{s - nL}{nL} \right),
         @f]
         where @f$a_n@f$ is the value of the function evaluated at the @f$ n @f$ knot and the other polynomial coefficients are determined by imposing continuity and 
         differentiability. 
         Parameter              | Interface parameter name             | Description 
         -----------------------|--------------------------------------|------------------------------------------------------------------------
         @f$\mathcal{R}(a_j)@f$ | <EM>particleName</EM>::Spline::Re::j | The real value of the lineshape evaluated at the @f$ j @f$th knot.
         @f$\mathcal{I}(a_j)@f$ | <EM>particleName</EM>::Spline::Im::j | The imaginary value of the lineshape evaluated at the @f$ j @f$th knot.
     */
    DECLARE_LINESHAPE( MIPWA );

    /** @ingroup Lineshapes class GSpline 
        @brief Lineshape with an arbitrary running width determined from a spline. 
      */
    DECLARE_LINESHAPE( GSpline );
    DECLARE_LINESHAPE( FormFactorSpline );
    DECLARE_LINESHAPE( DecaySpline );
    DECLARE_LINESHAPE( InelasticSpline );

    /** @ingroup Lineshapes class GenericKmatrix
        @brief Implementation of a generic K-matrix
      */
    DECLARE_LINESHAPE(GenericKmatrix);

    /** @ingroup Lineshapes class EtaDalitz 
        @brief Empirical Dalitz plot distribution for  @f$ \eta \rightarrow \pi^{+}\pi^{-}\pi^{0} @f$
        
        Lineshape that implements the Dalitz plot distribution for decays @f$ \eta \rightarrow \pi^{+}\pi^{-}\pi^{0} @f$, 
        parameterised in terms of the kinetic energy @f$T_x@f$ of each decay product in the rest frame of the @f$\eta@f$ meson.
        By convention, the third particle is neutral pion. The amplitude is given by: 
        @f[
          \mathcal{A}(T_1, T_2, T_3) = \sqrt{ 1 - 1.07\left( \frac{3 T_3}{T_1 + T_2 + T_3} -1 \right) } 
        @f]
        Multiplied by an arbitrary gaussian lineshape to account for the mass distribution in @f$s_{\pi^{+}\pi^{-}\pi^{0}}@f$.
      */   
    DECLARE_LINESHAPE( EtaDalitz );

    /** @ingroup Lineshapes class TD
        @brief (Linear) time dependence
        
        ''Lineshape'' that gives a linear time-dependence for use in time-dependent generation / fitting. 
        @f[
          \mathcal{A}(t) = \frac{t}{2\tau},
        @f]
        where @f$\tau@f$ is the proper decay time. 

      */
    DECLARE_LINESHAPE( TD );

    DECLARE_LINESHAPE( Photon );
  } // namespace Lineshape
  
  Expression Q2( const Expression& Msq, const Expression& M1sq, const Expression& M2sq );

  Expression kFactor( const Expression& mass, const Expression& width,
                      std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
  Expression BlattWeisskopf_Norm( const Expression& z, const Expression& z0, unsigned int L );

  Expression BL( const Expression& s, const Expression& s0, const Expression& s1, const Expression& s2,
                 const Expression& radius, const unsigned int& L );

  Expression BlattWeisskopf( const Expression& z, unsigned int L );

  Expression pol( const AmpGen::Expression& X, const std::vector<Expression>& p );

  std::vector<Expression> parameterVector(const std::string& name, const size_t& nParam);

  Expression width( const Expression& s, const Expression& s1, const Expression& s2, const Expression& mass,
                    const Expression& width, const Expression& radius, unsigned int L,
                    DebugSymbols* dbexpressions = nullptr );
} // namespace AmpGen
#endif
