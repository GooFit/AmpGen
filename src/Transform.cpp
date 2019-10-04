#include "AmpGen/Transform.h"

using namespace AmpGen; 

Transform::Transform( const Expression& arg, const Tensor& k, const Type& type ) :
  m_arg(arg),
  m_k(k),
  m_type(type) 
{}

Transform Transform::inverse() const 
{
  return Transform(m_arg,-m_k,m_type);
}

Expression AmpGen::cosh( const Expression& x ){ return (fcn::exp(x) + fcn::exp(-x))/2.;}
Expression AmpGen::sinh( const Expression& x ){ return (fcn::exp(x) - fcn::exp(-x))/2.;}
Expression AmpGen::tanh( const Expression& x ){ return (fcn::exp(x) - fcn::exp(-x))/(fcn::exp(x)+fcn::exp(-x));}

Expression AmpGen::acosh( const Expression& x ){ return fcn::log(x + fcn::sqrt(x*x-1));}
Expression AmpGen::asinh( const Expression& x ){ return fcn::log(x + fcn::sqrt(x*x+1));}
Expression AmpGen::atanh( const Expression& x ){ return 0.5*fcn::log((1+x)/(1-x));}

Tensor Transform::sigma_dot_p(const Tensor& p) const 
{
  return p[0] * Sigma[0] + p[1] * Sigma[1] + p[2] * Sigma[2];
}

Tensor Transform::J_dot_p( const Tensor& p ) const 
{
  Expression z(0);
  return Tensor(
      {     z, -p[2],  p[1] ,z 
      , p[2],     z, -p[0] ,z 
      ,-p[1],  p[0],     z ,z 
      ,    z,     z,     z ,z  }, Tensor::dim(4,4) );
}

Tensor Transform::K_dot_p( const Tensor& p ) const 
{
  Expression z(0);
  return Tensor(
      {    z,  z,  z ,p[0] 
      ,   z,  z,  z ,p[1] 
      ,   z,  z,  z ,p[2] 
      , p[0],  p[1],  p[2] ,z  }, Tensor::dim(4,4) );
}

Tensor Transform::make_bispinor_operator( const Tensor& A, const Tensor& B, const Tensor& C, const Tensor& D ) const
{
  return Tensor( {A[{0,0}], A[{0,1}], B[{0,0}],B[{0,1}],
      A[{1,0}], A[{1,1}], B[{1,0}],B[{1,1}],
      C[{0,0}], C[{0,1}], D[{0,0}],D[{0,1}],
      C[{1,0}], C[{1,1}], D[{1,0}],D[{1,1}] }, 
      Tensor::dim(4,4));
}

Tensor Transform::rotate_spinor() const 
{
  complex_t i (0,1);
  return Identity(2) * fcn::sqrt(0.5*(1+m_arg)) - i * sigma_dot_p(m_k) * fcn::sqrt(0.5*(1-m_arg));
}
Tensor Transform::boost_spinor() const 
{
  return Identity(2) * fcn::sqrt(0.5*(m_arg+1)) + sigma_dot_p(m_k)*fcn::sqrt(0.5*(m_arg-1));
}

Tensor Transform::operator()(const Representation& repr) const
{
  Tensor::Index m,j,k;
  Tensor I2 = Identity(2);
  Tensor I4 = Identity(4);
  if( m_type == Type::Rotate && repr == Representation::Spinor ) return rotate_spinor();
  if( m_type == Type::Boost  && repr == Representation::Spinor ) return boost_spinor();
  if( m_type == Type::Rotate && repr == Representation::Bispinor ){ 
    auto r = rotate_spinor();
    auto z = Tensor( Tensor::dim(2,2) );
    return make_bispinor_operator(r,z,z,r);
  }
  if( m_type == Type::Boost && repr == Representation::Bispinor )
  {
    Tensor t1 = I2 * fcn::sqrt(0.5*(m_arg+1)); 
    Tensor t2 = -sigma_dot_p(m_k) * fcn::sqrt(0.5*(m_arg-1));
    return make_bispinor_operator( t1, t2, t2, t1 );
  }
  if( m_type == Type::Rotate && repr == Representation::Vector )
  {
    auto K = J_dot_p(m_k); 
    return I4(j,k) + fcn::sqrt(1 - m_arg*m_arg)*K(j,k) + (1 - m_arg) * K(j,m) * K(m,k);
  }
  if( m_type == Type::Boost && repr == Representation::Vector )
  {
    auto K = K_dot_p(m_k); 
    return I4(j,k) - fcn::sqrt(m_arg*m_arg - 1)*K(j,k) + (m_arg - 1) * K(j,m) * K(m,k);
  }
  ERROR("Representation: " << repr << " " << m_type << " " << m_type << " not implemented");
  return Tensor();
}

TransformSequence TransformSequence::inverse() const
{
  TransformSequence rt; 
  for( auto i = m_transforms.rbegin(); i != m_transforms.rend(); ++i )
    rt.add( i->inverse() );
  return rt;
}

Tensor TransformSequence::operator()( const Transform::Representation& repr )
{
  if( m_cache[repr].nElements() != 1 )
  {
    return m_cache[repr];
  }
  if( m_transforms.size() == 0 ){
    if( repr == Transform::Representation::Spinor ) return Identity(2);
    else return Identity(4);
  }
  Tensor::Index a,b,c;
  Tensor rt = m_transforms[0](repr);
  rt.st();
  for( size_t i = 1 ; i < m_transforms.size(); ++i )
  {
    Tensor rti = m_transforms[i](repr);
    rti.st(true);
    rt = rti(a,b) * rt(b,c);
  }
  m_cache[repr] = rt;
  return rt;
}

Tensor TransformSequence::operator()( const Tensor& tensor, 
    const Transform::Representation& repr )
{
  Tensor::Index a,b,c;
  auto seq = this->operator()(repr);
  return seq(a,b)*tensor(b);
}

Tensor Transform::operator()( const Tensor& tensor, 
    const Transform::Representation& repr ) const
{
  Tensor::Index a,b;
  auto seq = this->operator()(repr);
  return seq(a,b)*tensor(b);
}

void TransformSequence::add( const Transform& transform )
{
  m_transforms.emplace_back( transform );
}

void TransformSequence::add( const TransformSequence& transform )
{
  for( auto& t : transform ) m_transforms.emplace_back(t);
}

void TransformSequence::clear()
{
  m_transforms.clear();
}

void TransformSequence::stepThrough(const Tensor& tensor, 
    const Transform::Representation& repr) 
{
  Tensor tr = tensor ;
  Tensor::Index a,b;
  INFO( "Frame[0] = " << tr[0]() << " " << tr[1]() << " " << tr[2]() << " " << tr[3]() );
  for(unsigned int i = 0 ; i < m_transforms.size(); ++i ){
    tr = m_transforms[i](repr)(a,b) * tr(b);
    INFO( "Frame["<<i+1<<"] = " << tr[0]() << " " << tr[1]() << " " << tr[2]() << " " << tr[3]() );
  }
}


