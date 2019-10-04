#include "AmpGen/Integrator.h"
#include <iostream>
using namespace AmpGen;
using namespace std;
Bilinears::Bilinears(const size_t& r, const size_t& c) : 
  rows(r), 
  cols(c), 
  norms(c*r), 
  markAsZero(c*r,0), 
  calculate(c*r,1) 
{}

complex_t Bilinears::get( const size_t& x, const size_t& y ) const { return norms[x * cols + y]; }
void      Bilinears::set( const size_t& x, const size_t& y,  const complex_t& f ) { norms[x * cols + y] = f; calculate[x*cols+y] = false; }
void Bilinears:: setZero( const size_t& x, const size_t& y ){ markAsZero[x * cols + y] = true; }
void Bilinears::resetCalculateFlags()
{ 
  for(size_t i = 0 ;i < calculate.size(); ++i) calculate[i] = true; 
} 
complex_t& Bilinears::operator()( const size_t& x, const size_t& y ){ return norms[x*cols+y];}
bool   Bilinears::isZero( const size_t& x, const size_t& y ){ return markAsZero[ x*cols+y] ; } 
bool Bilinears::workToDo( const size_t& x, const size_t& y ) const { 
  DEBUG("x = "<<x<<"\ny = "<<y<<"\ncols = "<<cols<<"\ncalculate size = "<<calculate.size());
  DEBUG("Calculate entry = "<<calculate[x * cols + y]);
  return calculate[x*cols+y] ; 
  } // && ! markAsZero[x*cols+y] ; }
void   Bilinears::resize( const size_t& r, const size_t& c)
{
  rows = r;
  cols = c;
  norms.resize( r * c );
  markAsZero.resize( r * c );
  calculate.resize(  r * c );
  for( size_t i = 0 ; i < r*c; ++i ){
    norms[i]      = 0;
    markAsZero[i] = false;
    calculate[i]  = true;
  }
}
