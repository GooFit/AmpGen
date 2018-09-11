#ifndef AMPGEN_LINESHAPES_H
#define AMPGEN_LINESHAPES_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"

namespace AmpGen
{
  struct ILineshape;
} // namespace AmpGen

#define DECLARE_LINESHAPE( X )                                                                                         \
  struct X : public AmpGen::ILineshape {                                                                               \
    static std::string _id;                                                                                            \
    X(){ DEBUG("Constructing lineshape") ;} \
    virtual AmpGen::Expression get( const AmpGen::Expression& s, const AmpGen::Expression& s1,                         \
                                    const AmpGen::Expression& s2, const std::string& particleName,                     \
                                    const unsigned int& L, const std::string& lineshapeModifier,                       \
                                    DebugSymbols* dbexpressions = 0 );                                                 \
  }

#define DECLARE_GENERIC_SHAPE( X )                                                                                     \
  struct X : public AmpGen::ILineshape {                                                                               \
    static std::string _id;                                                                                            \
    virtual AmpGen::Expression get( const std::vector<AmpGen::Tensor>& p, const std::string& particleName,             \
                                    const unsigned int& L, const std::string& lineshapeModifier,                       \
                                    AmpGen::DebugSymbols* dbexpressions = 0 );                                         \
  }

#define DEFINE_LINESHAPE( X )                                                                                          \
  REGISTER_WITH_KEY( ILineshape, Lineshape::X, #X, std::string );                                                      \
  AmpGen::Expression Lineshape::X::get( const AmpGen::Expression& s, const AmpGen::Expression& s1,                     \
                                        const AmpGen::Expression& s2, const std::string& particleName,                 \
                                        const unsigned int& L, const std::string& lineshapeModifier,                   \
                                        AmpGen::DebugSymbols* dbexpressions )

#define DEFINE_GENERIC_SHAPE( X )                                                                                      \
  REGISTER_WITH_KEY( ILineshape, Lineshape::X, #X, std::string );                                                      \
  AmpGen::Expression Lineshape::X::get( const std::vector<AmpGen::Tensor>& p, const std::string& particleName,         \
                                        const unsigned int& L, const std::string& lineshapeModifier,                   \
                                        AmpGen::DebugSymbols* dbexpressions )

namespace AmpGen
{
  class Tensor;

  /// some general helper functions ///
  Expression Q2( const Expression& Msq, const Expression& M1sq, const Expression& M2sq );

  Expression kFactor( const Expression& mass, const Expression& width,
                      std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
  Expression BlattWeisskopf_Norm( const Expression& z2, const Expression& z02, unsigned int L );

  Expression BL( const Expression& s, const Expression& s0, const Expression& s1, const Expression& s2,
                 const Expression& radius, const unsigned int& L );
  Expression BlattWeisskopf( const Expression& z2, unsigned int L );

  Expression pol( const AmpGen::Expression& X, const std::vector<Expression>& p );

  std::vector<Expression> parameterVector( const std::string& name, const unsigned int& nParam );

  Expression width( const Expression& s, const Expression& s1, const Expression& s2, const Expression& mass,
                    const Expression& width, const Expression& radius, unsigned int L,
                    DebugSymbols* dbexpressions = nullptr );

  /// Shape base class
  struct ILineshape {

    virtual ~ILineshape() = default;
    virtual Expression get( const Expression& s, const Expression& s1, const Expression& s2,
                            const std::string& particleName, const unsigned int& L,
                            const std::string& lineshapeModifier, DebugSymbols* dbexpressions = nullptr )
    {
      return 1;
    }

    virtual Expression get( const std::vector<AmpGen::Tensor>& p, const std::string& particleName,
                            const unsigned int& L, const std::string& lineshapeModifier,
                            AmpGen::DebugSymbols* dbexpressions = nullptr )
    {
      return 1;
    }
    ILineshape* create() { return this; }
  };

  /// static, global lineshape factory ///
  class LineshapeFactory : public Factory<ILineshape>
  {
  public:
    static Expression getLineshape( const std::string& lineshape, const Expression& s, const Expression& s1,
                                    const Expression& s2, const std::string& particleName, const unsigned int& L,
                                    std::vector<std::pair<std::string, Expression>>* dbexpressions = nullptr );
    static Expression getGenericShape( const std::string& lineshape,
                                       const std::vector<Tensor>& p, const std::string& particleName,
                                       const unsigned int& L, AmpGen::DebugSymbols* dbexpressions );
    static bool isLineshape( const std::string& lineshape );
  };

  namespace Lineshape
  {
    /// Two-body Breit-Wigner lineshape
    DECLARE_LINESHAPE( BW );
    /// Breit-Wigner lineshape with fixed width
    DECLARE_LINESHAPE( SBW );
    /// Non-relativistic Breit-Wigner lineshape
    DECLARE_LINESHAPE( NonRelBW );

    /// Gounaris-Sakurai lineshape models dispersive corrections to the I=1,J=1 \f$\pi\pi\f$ propagator (G.J. Gounaris and J.J. Sakurai, Phys. Rev. Lett.21, 244 (1968))
    DECLARE_LINESHAPE( GounarisSakurai );

    /// LASS shape used to model the \f$ K\pi \f$ S-wave
    DECLARE_LINESHAPE( LASS );

    DECLARE_LINESHAPE( gLASS );
    /// Flatte lineshape to describe resonances with coupled channels such as \f$f_{0}(980)^{0} / a_{0}(980) \f$ (S.M.Flatté, Phys. Lett B. 63, 224 (1976))
    DECLARE_LINESHAPE( Flatte );
    DECLARE_LINESHAPE( Bugg );
    DECLARE_LINESHAPE( Isotensor );

    /// Exponential form-factor of type \f$ e^{ - q^2 R^2 } \f$
    DECLARE_LINESHAPE( ExpFF );

    DECLARE_LINESHAPE( FormFactor );
    DECLARE_LINESHAPE( NonRes );
    DECLARE_LINESHAPE( None );

    /// Gaussian lineshape \f$ = e^{ -(x-\mu)^2 / 2\sigma^{2} } \f$
    DECLARE_LINESHAPE( Gaussian );

    /// Polynominal shape
    DECLARE_LINESHAPE( Poly );

    /// Anisovich-Sarantsev Isoscalar K-matrix from https://arxiv.org/abs/hep-ph/0204328
    DECLARE_LINESHAPE( kMatrix );

    /// I=1/2 and I=3/2 K Matrices used in the description of the \f$ K\pi \f$ S-wave in the analysis of \f$ D^{+}
    /// \rightarrow K^{-}\pi^{+}\pi^{+}\f$ https://arxiv.org/abs/0705.2248
    DECLARE_LINESHAPE( FOCUS );

    /// I=1/2 and I=3/2 K Matrices used in the description of the \f$ K\pi\f$ S-wave in the analysis of \f$\eta_{c}
    /// \rightarrow K^{+} K^{-} \pi^{0}\f$ https://arxiv.org/abs/1701.04881
    DECLARE_LINESHAPE( PALANO );

    /// I=1, J=1 K Matrix used to describe the \f$ \rho(770), \rho(1450), \rho(1900) \f$ system, including the \f$
    /// \pi\pi,  KK , \pi\pi\pi\pi \f$ channels
    DECLARE_LINESHAPE( ObelixRho );

    /// K matrix to describe \f$K_1(1270) / K_1(1400)\f$, WARNING, does not work at intended. 
    DECLARE_LINESHAPE( AxialKaon );
    /// Flexible K matrix implementation that accepts a variable number of scalar channels and pole terms
    DECLARE_LINESHAPE( kMatrixSimple );

    /// Spline based lineshapes
    DECLARE_LINESHAPE( MIPWA );
    DECLARE_LINESHAPE( WidthSpline );
    DECLARE_LINESHAPE( GSpline );
    DECLARE_LINESHAPE( FormFactorSpline );
    DECLARE_LINESHAPE( DecaySpline );
    DECLARE_LINESHAPE( InelasticSpline );

    /// Implements Dalitz plot distribution for decays \f$ \eta \rightarrow \pi^{+}\pi^{-}\pi^{0}\f$
    DECLARE_GENERIC_SHAPE( EtaDalitz );

  } // namespace Lineshape
} // namespace AmpGen
#endif
