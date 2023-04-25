#include "AmpGen/Vertex.h"
#include "AmpGen/DiracMatrices.h"
#include "AmpGen/Tensor.h"
using namespace AmpGen; 

namespace { 
  static const Tensor::Index mu    = Tensor::Index();
  static const Tensor::Index nu    = Tensor::Index();
  static const Tensor::Index alpha = Tensor::Index();
  static const Tensor::Index beta  = Tensor::Index();
  static const Tensor::Index a     = Tensor::Index();
  static const Tensor::Index b     = Tensor::Index(); 
  static const Tensor::Index c     = Tensor::Index(); 
  static const Tensor::Index d     = Tensor::Index(); 
}

DEFINE_VERTEX( f_fS_S )
{
  Tensor f_fS_Sv = Spin1hProjector(P)(a,b) * V1(b) * V2[0];
  ADD_DEBUG_TENSOR( f_fS_Sv, db);
  return f_fS_Sv; 
}

DEFINE_VERTEX( f_fS_P )
{
  Tensor proj   = Spin1hProjector(P);
  Tensor L      = Orbital_PWave(P,Q);
  ADD_DEBUG_TENSOR( P, db );
  ADD_DEBUG_TENSOR( Q, db );
  ADD_DEBUG_TENSOR( L, db );
  Tensor t      = proj(a, b) * Gamma[4](b,c) * slash(L)(c,d) * V1(d);
  return t; 
}

DEFINE_VERTEX( f_Vf_S )
{
  Tensor proj   = Spin1hProjector(P);
  return proj(a, b) * Gamma[4](b,c) * gamma_twiddle(P)(mu,c,d) * V2(d) * V1(-mu);
}

DEFINE_VERTEX( f_Vf_P ) /// = A3 
{
  Tensor proj   = Spin1hProjector(P);
  auto L = Orbital_PWave(P,Q);
  return proj( alpha, beta ) * V2(beta) * dot( V1, L );
}

DEFINE_VERTEX( f_Vf_P1 ) //// = A6
{
  Tensor proj   = Spin1hProjector(P);
  Tensor L = Orbital_PWave(P,Q);
  Tensor t = LeviCivita()(mu,nu,alpha,beta) * V1(-nu) * L(-alpha) * P(-beta); 
  t.st();
  return proj( a, b ) * Gamma[4](b,c) * gamma_twiddle(P)(mu,c,d) * V2(d) * t(-mu);
}

DEFINE_VERTEX( f_Vf_D ) //// = A8 
{
  Tensor proj   = Spin1hProjector(P);
  Tensor gt     = gamma_twiddle(P);
  Tensor L      = Orbital_DWave(P,Q)(-mu,-nu) * V1(nu); 
  return proj( a, b ) * Gamma[4](b,c) * gt(mu,c,d) * V2(d) * L(-mu);
}

DEFINE_VERTEX( f_Vf_S1 )
{
  Tensor proj   = Spin1hProjector(P);
  Tensor vSlash = gamma_twiddle(P)(mu,a,b) * V1(-mu);
  return proj( a, b ) * vSlash(b,c) * V2(c);
}

DEFINE_VERTEX( f_Vf_SL )
{ 
  Tensor proj = Spin1hProjector(P);
  return proj(a, b) * Gamma4Vec()(mu,b,c) * ( Identity(4) - Gamma[4] )(c,d)* V2(d) * V1(-mu);
}

DEFINE_VERTEX( f_Vf_SR )
{ 
  Tensor proj = Spin1hProjector(P);
  return proj(a, b) * Gamma4Vec()(mu,b,c) * ( Identity(4) + Gamma[4] )(c,d)* V2(d) * V1(-mu);
}

DEFINE_VERTEX( f_fS_SL )
{ 
  Tensor proj = Spin1hProjector(P);
  return proj(a, b) * ( Identity(4) - Gamma[4] )(b,c)* V2(c);
}


DEFINE_VERTEX( f_fS_SR )
{ 
  Tensor proj = Spin1hProjector(P);
  return proj(a, b) * ( Identity(4) + Gamma[4] )(b,c)* V2(c);
}

DEFINE_VERTEX( f_fS_S1 ){ return Spin1hProjector(P)(a,b) * Gamma[4](b,c) * V1(c) * V2[0];}

DEFINE_VERTEX( f_fS_P1 )
{
  return Spin1hProjector(P)(a, b) * slash(Orbital_PWave(P,Q))(b,c) * V1(c) * V2[0];
}

DEFINE_VERTEX( f_Tf_P )
{
  Tensor proj   = Spin1hProjector(P);
  Tensor T = V1;
  T.imposeSymmetry(0,1);
  T.st();
  return proj( a, b ) * gamma_twiddle(P)(mu,b,c) * V2(c) * T(-mu,-nu) * Orbital_PWave(P,Q)(nu);
}

DEFINE_VERTEX( f_Vf_P2 ) //A4 
{
  Tensor proj   = Spin1hProjector(P);
  return proj( alpha, beta ) * Gamma[4](beta,nu) * V2(nu) * dot( V1, Orbital_PWave( P, Q ) );
}

DEFINE_VERTEX( f_Vf_P3 ) //A5
{
  Tensor proj   = Spin1hProjector(P);
  Tensor L = Orbital_PWave(P,Q);
  Tensor t = LeviCivita()(-mu,-nu,-alpha,-beta) * L(nu) * V1(alpha) * P(beta); 
  t.st();
  return proj( a, b ) * gamma_twiddle(P)(mu,b,c) * V2(c) * t(-mu); 
}


DEFINE_VERTEX( f_Vf_D1 )
{
  Tensor proj   = Spin1hProjector(P);
  Tensor gt     = gamma_twiddle(P);
  Tensor L      = Orbital_DWave(P,Q)(-mu,-nu) * V1(nu); 
  return proj( a, b ) * gt(mu,b,c) * L(-mu) * V2(c) ;
}


DEFINE_VERTEX( r_fS_P )
{
  Tensor L = Orbital_PWave(P,Q);
  L.st();
  Tensor F = Spin1hProjector(P);
  Tensor V = ( (-1) * L(nu) * F(a,d) + (1./3.) * F(a,b) * gamma_twiddle(P)(nu,b,c) * slash(L)(c,d) ) * V1(d); 
  V.st();
  return V;
}

DEFINE_VERTEX( r_fS_D )
{
  Tensor::Index e;
  Tensor L = Orbital_PWave(P,Q);
  Tensor F = Spin1hProjector(P);
  L.st(true);
  Expression L2 = make_cse(dot(L,L));
  Tensor gt = gamma_twiddle(P);
  Tensor rt =  ( L(mu) * F(a,b) * slash(L)(b,c) - (L2/3.) * F(a,b) * gt(mu,b,c) ) * Gamma[4](c,d) * V1(d); 
  rt.st();
  return rt;
}


DEFINE_VERTEX( f_rS_D )
{
  Tensor F = Spin1hProjector(P)(a,b) 
    * Gamma[4](b,d) 
    * gamma_twiddle(P)(mu,d,c)
    * V1(alpha,c) 
    * Orbital_DWave(P,Q)(-mu,-alpha) 
    * V2[0];
  F.st();
  return F;
}

DEFINE_VERTEX( f_rS_P )
{
  auto L = Orbital_PWave(P,Q);
  Tensor F = Spin1hProjector(P)(a,b)                 * V1(-mu,b) * L(mu) * V2[0];
  ADD_DEBUG_TENSOR( V1, db );
  ADD_DEBUG_TENSOR( L , db );
  ADD_DEBUG_TENSOR( F , db );
  F.st();
  return F;
}

DEFINE_VERTEX( f_rS_P1 )
{
  auto L = Orbital_PWave(P,Q);
  Tensor F = Spin1hProjector(P)(a,b) * Gamma[4](b,c) * V1(-mu,c) * L(mu) * V2[0]; 
  F.st();
  return F;
}

DEFINE_VERTEX( F_FS_S )
{
  Tensor f_fS_Sv = Spin1hbProjector(P)(a,b) * V1(b) * V2[0];
  ADD_DEBUG_TENSOR( f_fS_Sv, db);
  return f_fS_Sv; 
}

DEFINE_VERTEX( F_FS_P )
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor L      = Orbital_PWave(P,Q);
  ADD_DEBUG_TENSOR( P, db );
  ADD_DEBUG_TENSOR( Q, db );
  ADD_DEBUG_TENSOR( L, db );
  Tensor t      = proj(a, b) * Gamma[4](b,c) * slash(L)(c,d) * V1(d);
  return t; 
}

DEFINE_VERTEX( F_VF_S )
{
  Tensor proj   = Spin1hbProjector(P);
  return proj(a, b) * Gamma[4](b,c) * gamma_twiddle(P)(mu,c,d) * V2(d) * V1(-mu);
}

DEFINE_VERTEX( F_VF_P ) /// = A3 
{
  Tensor proj   = Spin1hbProjector(P);
  auto L = Orbital_PWave(P,Q);
  return proj( alpha, beta ) * V2(beta) * dot( V1, L );
}

DEFINE_VERTEX( F_VF_P1 ) //// = A6
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor L = Orbital_PWave(P,Q);
  Tensor t = LeviCivita()(mu,nu,alpha,beta) * V1(-nu) * L(-alpha) * P(-beta); 
  t.st();
  return proj( a, b ) * Gamma[4](b,c) * gamma_twiddle(P)(mu,c,d) * V2(d) * t(-mu);
}

DEFINE_VERTEX( F_VF_D ) //// = A8 
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor gt     = gamma_twiddle(P);
  Tensor L      = Orbital_DWave(P,Q)(-mu,-nu) * V1(nu); 
  return proj( a, b ) * Gamma[4](b,c) * gt(mu,c,d) * V2(d) * L(-mu);
}

DEFINE_VERTEX( F_VF_S1 )
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor vSlash = gamma_twiddle(P)(mu,a,b) * V1(-mu);
  return proj( a, b ) * vSlash(b,c) * V2(c);
}

DEFINE_VERTEX( F_VF_SL )
{ 
  Tensor proj = Spin1hbProjector(P);
  return proj(a, b) * Gamma4Vec()(mu,b,c) * ( Identity(4) - Gamma[4] )(c,d)* V2(d) * V1(-mu);
}

DEFINE_VERTEX( F_VF_SR )
{ 
  Tensor proj = Spin1hbProjector(P);
  return proj(a, b) * Gamma4Vec()(mu,b,c) * ( Identity(4) + Gamma[4] )(c,d)* V2(d) * V1(-mu);
}

DEFINE_VERTEX( F_FS_SL )
{ 
  Tensor proj = Spin1hbProjector(P);
  return proj(a, b) * ( Identity(4) - Gamma[4] )(b,c)* V2(c);
}


DEFINE_VERTEX( F_FS_SR )
{ 
  Tensor proj = Spin1hbProjector(P);
  return proj(a, b) * ( Identity(4) + Gamma[4] )(b,c)* V2(c);
}

DEFINE_VERTEX( F_FS_S1 ){ return Spin1hProjector(P)(a,b) * Gamma[4](b,c) * V1(c) * V2[0];}

DEFINE_VERTEX( F_FS_P1 )
{
  return Spin1hbProjector(P)(a, b) * slash(Orbital_PWave(P,Q))(b,c) * V1(c) * V2[0];
}

DEFINE_VERTEX( F_TF_P )
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor T = V1;
  T.imposeSymmetry(0,1);
  T.st();
  return proj( a, b ) * gamma_twiddle(P)(mu,b,c) * V2(c) * T(-mu,-nu) * Orbital_PWave(P,Q)(nu);
}

DEFINE_VERTEX( F_VF_P2 ) //A4 
{
  Tensor proj   = Spin1hbProjector(P);
  return proj( alpha, beta ) * Gamma[4](beta,nu) * V2(nu) * dot( V1, Orbital_PWave( P, Q ) );
}

DEFINE_VERTEX( F_VF_P3 ) //A5
{
  Tensor proj   = Spin1hbProjector(P);
  Tensor L = Orbital_PWave(P,Q);
  Tensor t = LeviCivita()(-mu,-nu,-alpha,-beta) * L(nu) * V1(alpha) * P(beta); 
  t.st();
  return proj( a, b ) * gamma_twiddle(P)(mu,b,c) * V2(c) * t(-mu); 
}


DEFINE_VERTEX( F_VF_D1 )
{
  Tensor proj   = Spin1hProjector(P);
  Tensor gt     = gamma_twiddle(P);
  Tensor L      = Orbital_DWave(P,Q)(-mu,-nu) * V1(nu); 
  return proj( a, b ) * gt(mu,b,c) * L(-mu) * V2(c) ;
}


DEFINE_VERTEX( R_FS_P )
{
  Tensor L = Orbital_PWave(P,Q);
  L.st();
  Tensor F = Spin1hbProjector(P);
  Tensor V = ( (-1) * L(nu) * F(a,d) + (1./3.) * F(a,b) * gamma_twiddle(P)(nu,b,c) * slash(L)(c,d) ) * V1(d); 
  V.st();
  return V;
}

DEFINE_VERTEX( R_FS_D )
{
  Tensor::Index e;
  Tensor L = Orbital_PWave(P,Q);
  Tensor F = Spin1hbProjector(P);
  L.st(true);
  Expression L2 = make_cse(dot(L,L));
  Tensor gt = gamma_twiddle(P);
  Tensor rt =  ( L(mu) * F(a,b) * slash(L)(b,c) - (L2/3.) * F(a,b) * gt(mu,b,c) ) * Gamma[4](c,d) * V1(d); 
  rt.st();
  return rt;
}


DEFINE_VERTEX( F_RS_D )
{
  Tensor F = Spin1hbProjector(P)(a,b) 
    * Gamma[4](b,d) 
    * gamma_twiddle(P)(mu,d,c)
    * V1(alpha,c) 
    * Orbital_DWave(P,Q)(-mu,-alpha) 
    * V2[0];
  F.st();
  return F;
}

DEFINE_VERTEX( F_RS_P )
{
  auto L = Orbital_PWave(P,Q);
  Tensor F = Spin1hbProjector(P)(a,b)                 * V1(-mu,b) * L(mu) * V2[0];
  ADD_DEBUG_TENSOR( V1, db );
  ADD_DEBUG_TENSOR( L , db );
  ADD_DEBUG_TENSOR( F , db );
  F.st();
  return F;
}

DEFINE_VERTEX( F_RS_P1 )
{
  auto L = Orbital_PWave(P,Q);
  Tensor F = Spin1hbProjector(P)(a,b) * Gamma[4](b,c) * V1(-mu,c) * L(mu) * V2[0]; 
  F.st();
  return F;
}
