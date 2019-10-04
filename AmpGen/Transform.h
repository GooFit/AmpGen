#ifndef AMPGEN_TRANSFORM_H
#define AMPGEN_TRANSFORM_H

#include "AmpGen/Expression.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/DiracMatrices.h"

namespace AmpGen {
  Expression cosh(const Expression& x);
  Expression sinh(const Expression& x);
  Expression tanh(const Expression& x);

  Expression acosh(const Expression& x); 
  Expression asinh(const Expression& x); 
  Expression atanh(const Expression& x); 

  class Transform 
  {
    public:
      enum Type           { Boost, Rotate };
      enum Representation { Spinor, Bispinor, Vector }; 

      Transform( const Expression& arg, const Tensor& k, const Type& type );
      Transform inverse() const; 
      Tensor operator()(const Representation& repr) const;
      Tensor operator()(const Tensor& tensor, const Representation& repr=Representation::Vector) const; 
    private:

      Tensor sigma_dot_p(const Tensor& p) const;
      Tensor J_dot_p( const Tensor& p ) const;
      Tensor K_dot_p( const Tensor& p ) const;
      Tensor make_bispinor_operator( const Tensor& A, const Tensor& B, const Tensor& C, const Tensor& D ) const;
      Tensor rotate_spinor() const;
      Tensor boost_spinor() const;

      Expression m_arg; 
      Tensor     m_k;
      Type       m_type; 
  };

  class TransformSequence
  {
    public:
      TransformSequence() = default; 
      TransformSequence inverse() const;
      Tensor operator()( const Transform::Representation& repr );
      Tensor operator()( const Tensor& tensor, 
          const Transform::Representation& repr=Transform::Representation::Vector );
      void add( const Transform& transform );
      void add( const TransformSequence& transform );
      void stepThrough( const Tensor& tensor, 
                        const Transform::Representation& repr = Transform::Representation::Vector );
      
      void clear(); 
      std::vector<Transform>::reverse_iterator rbegin()       { return m_transforms.rbegin(); }
      std::vector<Transform>::reverse_iterator rend()         { return m_transforms.rend(); }
      std::vector<Transform>::const_iterator   begin()  const { return m_transforms.cbegin(); }
      std::vector<Transform>::const_iterator   end()    const { return m_transforms.cend(); }
      std::vector<Transform>::iterator         begin()        { return m_transforms.begin(); }
      std::vector<Transform>::iterator         end()          { return m_transforms.end(); }
    private:
      std::vector<Transform> m_transforms;
      std::array<Tensor,3>   m_cache; 
  };

}

#endif
