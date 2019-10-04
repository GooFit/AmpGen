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
  if ( Y != nullptr ) for( unsigned i = 0 ; i < X.size(); ++i ) \
    Y->emplace_back( std::string(#X) + Tensor::coordinates_to_string( X.coords(i) ) , X[i] );

#define ADD_DEBUG_TENSOR_NAMED( X, Y, Z )  \
  if ( Y != nullptr ) for( unsigned i = 0 ; i < X.size(); ++i ) \
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
      public:
        bool operator==( const Tensor::Index& other ) const { return m_ptr.get() == other.m_ptr.get(); }
        bool operator!=( const Tensor::Index& other ) const { return m_ptr.get() != other.m_ptr.get(); }
        bool isUpper() const { return m_isUpper; };
        explicit Index( bool isUpper = false ) : m_ptr( new int() ), m_isUpper( isUpper ) {}
        Index( const std::shared_ptr<int>& index, bool isUpper = false ) : m_ptr( index ), m_isUpper( isUpper ) {}
        Index operator-() const { return Index( m_ptr, !m_isUpper ); }
        friend std::ostream& operator <<( std::ostream& out, const Index& index );  
    };

    Tensor();
    explicit Tensor(const std::vector<Expression>& elements);
    explicit Tensor(const std::vector<unsigned>& dim);

    template <class TYPE> Tensor(const std::initializer_list<TYPE>& elements, 
                                 const std::vector<unsigned>& dim) : m_dim(dim)
    {
      setupCoordinates();
      for ( auto& x : elements ) append( x );
    }

    template <class TYPE> Tensor(const std::vector<TYPE>& elements, 
                                 const std::vector<unsigned>& dim) : m_dim(dim)
    {
      setupCoordinates();
      for ( auto& x : elements ) append( x );
    }

    /// Low level access of elements, either by coordinates or by index ///
    Expression& operator[]( const unsigned& i ); 
    Expression& operator[]( const std::vector<unsigned>& co ); 
    const Expression& operator[]( const unsigned& i ) const; 
    const Expression& operator[]( const std::vector<unsigned>& co ) const; 

    Expression get( const unsigned& co );
    Expression get( const unsigned& co ) const;
    Expression get( const std::vector<unsigned>& _co ) const;

    /// TensorProxy access to class members
    /// High level access is done via these commands, i.e. () operators

    Expression& operator()( const unsigned& a ) { return Tensor::operator[]( {a} ) ; }
    Expression& operator()( const unsigned& a, const unsigned& b) { return Tensor::operator[]( {a,b} ) ; }
    const Expression  operator()( const unsigned& a ) const { return Tensor::operator[]( {a} ) ; }
    const Expression& operator()( const unsigned& a, const unsigned& b) const { return Tensor::operator[]( {a,b} ) ; }

    TensorProxy operator()( const Tensor::Index& a ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c ) const;
    TensorProxy operator()( const Tensor::Index& a, const Tensor::Index& b, const Tensor::Index& c,
                             const Tensor::Index& d ) const;
    TensorProxy operator()( const std::vector<Tensor::Index>& indices ) const;
    
    Tensor operator-() const; 

    void st(const bool simplify=false);
    bool rankMatches( const Tensor& other );
   
    void imposeSymmetry( unsigned indexA, unsigned indexB);
    void imposeSymmetry( std::vector<unsigned> indices );

    Tensor Invert() const;
    std::string to_string(const ASTResolver* resolver=nullptr) const ;

    int metricSgn( const std::vector<unsigned>& coordinates ) const;
    int metricSgn( const unsigned& index ) const;
    void append( const Expression& expression );
    void append( const real_t& value );
    void append( const complex_t& value );
    void append( const std::string& value );
    void setupCoordinates(); 

    unsigned nDim() const;
    unsigned rank() const;
    unsigned size() const;
    unsigned index( const std::vector<unsigned>& _co ) const;
    unsigned symmetrisedIndex( const std::vector<unsigned>& _co ) const;
    unsigned nElements() const;

    const std::vector<unsigned> coords( const unsigned& index ) const;
    const std::vector<unsigned>& dims() const { return m_dim; }
    const std::string dimString() const ;

    void print(const bool& eval = false) const;

    const std::vector<unsigned>& uniqueElements() const { return m_uniqueElements; }
    void operator+=( const Tensor& rhs );
    void operator-=( const Tensor& rhs );
    Tensor conjugate() const;

    static std::vector<unsigned> index_to_coordinates( const unsigned& index, const std::vector<unsigned>& dim );
    static unsigned coordinates_to_index( const std::vector<unsigned>& coords, const std::vector<unsigned>& dim );
    static std::string coordinates_to_string( const std::vector<unsigned >& coordinates );
    template <class... ARGS>
    static std::vector<unsigned> dim( const ARGS&... args ){
      std::vector<unsigned> rt; 
      auto up = std::tuple<ARGS...>(args...);
      for_each(up, [&rt]( const unsigned& f ) { rt.emplace_back(f); } );
      return rt;
    }
  private:
    std::vector<unsigned> m_dim;   
    std::vector<unsigned> m_symmetrisedCoordinates; 
    std::vector<unsigned> m_uniqueElements; 
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

  Tensor operator+(const Tensor&, const Tensor&);
  Tensor operator-(const Tensor&, const Tensor&);
  Tensor operator/(const Tensor&, const Expression&);
  Tensor operator*(const Expression&, const Tensor&);
  Tensor operator*(const Tensor&, const Expression&);

  Tensor operator/(const Tensor&, const double&);
  Tensor operator*(const double&, const Tensor&);
  Tensor operator*(const Tensor&, const double&);

  TensorProxy operator*(const TensorProxy&, const TensorProxy&);
  TensorProxy operator+(const TensorProxy&, const TensorProxy&);
  TensorProxy operator-(const TensorProxy&, const TensorProxy&);

  TensorProxy operator/(const TensorProxy&, const Expression& );
  TensorProxy operator*(const Expression& , const TensorProxy&);
  TensorProxy operator*(const TensorProxy&, const Expression& );

  TensorProxy operator/(const TensorProxy&, const double&);
  TensorProxy operator*(const double&     , const TensorProxy&);
  TensorProxy operator*(const TensorProxy&, const double&);

  Tensor Identity( const unsigned& rank = 4 );
  
  const Tensor LeviCivita( const unsigned& rank = 4);
  Expression dot( const Tensor& A, const Tensor& B );

  std::ostream& operator <<( std::ostream& out, const Tensor::Index& index );
} // namespace AmpGen

#endif
