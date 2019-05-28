#ifndef AMPGEN_TENSOR_H
#define AMPGEN_TENSOR_H
#include <memory.h>
#include <stddef.h>
#include <algorithm>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "AmpGen/Expression.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Types.h"

#define ADD_DEBUG_TENSOR( X, Y )  \
  if ( Y != nullptr ) for( size_t i = 0 ; i < X.size(); ++i ) \
    Y->emplace_back( std::string(#X) + Tensor::coordinates_to_string( X.coords(i) ) , X[i] );

#define ADD_DEBUG_TENSOR_NAMED( X, Y, Z )  \
  if ( Y != nullptr ) for( size_t i = 0 ; i < X.size(); ++i ) \
    Y->emplace_back( Z + Tensor::coordinates_to_string( X.coords(i) ) , X[i] );

namespace AmpGen
{
  class ASTResolver;
  class TensorProxy;
  
  class Tensor
  {
  public:
    class Index {
      private: 
        std::shared_ptr<int> m_ptr;
        bool m_isUpper;
        std::string name;
      public:
        bool operator==( const Tensor::Index& other ) const { return m_ptr.get() == other.m_ptr.get(); }
        bool operator!=( const Tensor::Index& other ) const { return m_ptr.get() != other.m_ptr.get(); }
        bool isUpper() const { return m_isUpper; };
        Index( bool isUpper = false ) : m_ptr( new int() ), m_isUpper( isUpper ) {}
        Index( const std::shared_ptr<int>& index, bool isUpper = false ) : m_ptr( index ), m_isUpper( isUpper ) {}
        Index operator-() const { return Index( m_ptr, !m_isUpper ); }
        friend std::ostream& operator <<( std::ostream& out, const Index& index );  
    };

    Tensor();
    Tensor(const std::vector<Expression>& elements);
    Tensor(const std::vector<size_t>& dim);

    template <class TYPE> Tensor(const std::initializer_list<TYPE>& elements, 
                                 const std::vector<size_t>& dim) : m_dim(dim)
    {
      setupCoordinates();
      for ( auto& x : elements ) append( x );
    }

    template <class TYPE> Tensor(const std::vector<TYPE>& elements, 
                                 const std::vector<size_t>& dim) : m_dim(dim)
    {
      setupCoordinates();
      for ( auto& x : elements ) append( x );
    }

    /// Low level access of elements, either by coordinates or by index ///
    Expression& operator[]( const size_t& i ); 
    Expression& operator[]( const std::vector<size_t>& co ); 
    const Expression& operator[]( const size_t& i ) const; 
    const Expression& operator[]( const std::vector<size_t>& co ) const; 

    Expression get( const size_t& co );
    Expression get( const size_t& co ) const;
    Expression get( const std::vector<size_t>& _co ) const;

    /// TensorProxy access to class members
    /// High level access is done via these commands, i.e. () operators

    Expression& operator()( const size_t& a ) { return Tensor::operator[]( {a} ) ; }
    Expression& operator()( const size_t& a, const size_t& b) { return Tensor::operator[]( {a,b} ) ; }
    const Expression  operator()( const size_t& a ) const { return Tensor::operator[]( {a} ) ; }
    const Expression& operator()( const size_t& a, const size_t& b) const { return Tensor::operator[]( {a,b} ) ; }

    TensorProxy operator()( const Tensor::Index& a ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c,
                             const Tensor::Index& d ) const;
    TensorProxy operator()( const std::vector<Tensor::Index>& indices ) const;
    
    Tensor operator-() const; 

    void st(const bool simplify=false);
    bool rankMatches( const Tensor& other );
   
    void imposeSymmetry( size_t indexA, size_t indexB);
    void imposeSymmetry( std::vector<size_t> indices );

    Tensor Invert() const;
    std::string to_string(const ASTResolver* resolver=nullptr) const ;

    int metricSgn( const std::vector<size_t>& coordinates ) const;
    int metricSgn( const size_t& index ) const;
    void append( const Expression& expression );
    void append( const real_t& value );
    void append( const complex_t& value );
    void append( const std::string& value );
    void setupCoordinates(); 

    size_t nDim() const;
    size_t rank() const;
    size_t size() const;
    size_t index( const std::vector<size_t>& _co ) const;
    size_t symmetrisedIndex( const std::vector<size_t>& _co ) const;
    size_t nElements() const;

    const std::vector<size_t> coords( const size_t& index ) const;
    const std::vector<size_t>& dims() const { return m_dim; }
    const std::string dimString() const ;

    void print(const bool& eval = false) const;

    const std::vector<size_t>& uniqueElements() const { return m_uniqueElements; }
    void operator+=( const Tensor& rhs );
    void operator-=( const Tensor& rhs );
    Tensor conjugate() const;

    static std::vector<size_t> index_to_coordinates( const size_t& index, const std::vector<size_t>& dim );
    static size_t coordinates_to_index( const std::vector<size_t>& coords, const std::vector<size_t>& dim );
    static std::string coordinates_to_string( const std::vector<size_t >& coordinates );
    template <class... ARGS>
    static std::vector<size_t> dim( const ARGS&... args ){
      std::vector<size_t> rt; 
      auto up = std::tuple<ARGS...>(args...);
      for_each(up, [&rt]( const size_t& f ) { rt.emplace_back(f); } );
      return rt;
    }
  private:
    std::vector<size_t> m_dim;   
    std::vector<size_t> m_symmetrisedCoordinates; 
    std::vector<size_t> m_uniqueElements; 
    std::vector<Expression> m_elements;
  };

  /** @class TensorProxy
      Utility class that wraps a tensor and a set of indices such
      that tensor operations can be performed.
      It normally isn't useful to use this class directly,
      instead is designed as an intermediate proxy object for performing
      tensor manipulations
   */
  class TensorProxy
  {
  public:
    TensorProxy( const Tensor& tensor, const std::vector<Tensor::Index>& indices ) ;
    std::vector<Tensor::Index> indices() const;
    const Tensor& tensor() const;
    Tensor& tensorMutable();
    TensorProxy reorder( const std::vector<Tensor::Index>& indices );
    operator Expression() const;
    operator Tensor() const;
  
  private:
    Tensor m_tensor;
    std::vector<Tensor::Index> m_indices;
  };
  
  class TensorExpression : public IExpression 
  {
    public:
      TensorExpression( const Tensor& tensor );
      std::string to_string(const ASTResolver* resolver) const override;
      void resolve( ASTResolver& resolver ) const override; 
      complex_t operator()() const override;
      operator Expression() const;
    
    private:
      Tensor m_tensor; 
  };

  Tensor operator+( const Tensor& t1, const Tensor& t2 );
  Tensor operator-( const Tensor& t1, const Tensor& t2 );
  Tensor operator/( const Tensor& t1, const Expression& t2 );
  Tensor operator*( const Expression& t1, const Tensor& t2 );
  Tensor operator*( const Tensor& t1, const Expression& t2 );

  Tensor operator/( const Tensor& t1, const double& t2 );
  Tensor operator*( const double& t1, const Tensor& t2 );
  Tensor operator*( const Tensor& t1, const double& t2 );

  TensorProxy operator*( const TensorProxy& t1, const TensorProxy& t2 );
  TensorProxy operator+( const TensorProxy& t1, const TensorProxy& t2 );
  TensorProxy operator-( const TensorProxy& t1, const TensorProxy& t2 );

  TensorProxy operator/( const TensorProxy& t1, const Expression& t2 );
  TensorProxy operator*( const Expression& t1, const TensorProxy& t2 );
  TensorProxy operator*( const TensorProxy& t1, const Expression& t2 );

  TensorProxy operator/( const TensorProxy& t1, const double& t2 );
  TensorProxy operator*( const double& t1, const TensorProxy& t2 );
  TensorProxy operator*( const TensorProxy& t1, const double& t2 );

  Tensor Identity( const size_t& rank = 4 );
  
  const Tensor LeviCivita( const size_t& rank = 4);
  Expression dot( const Tensor& A, const Tensor& B );

  std::ostream& operator <<( std::ostream& out, const Tensor::Index& index );
} // namespace AmpGen

#endif
