
#include <gsl/gsl_integration.h>
#include <iostream>

#include "AmpGen/MetaUtils.h"
#include "AmpGen/simd/utils.h"
#include "AmpGen/simd/integrate_fp.h"
#include <queue>

namespace AmpGen {
  template <typename function_type> double integrate_1d( const function_type& fcn, const double& min, const double& max, gsl_integration_workspace* ws = nullptr)
  {
    gsl_integration_workspace* w = (ws == nullptr) ? gsl_integration_workspace_alloc (1000) : ws;
    double result = 0;
    double error  = 0;
    gsl_function F;
    //  static_assert( is_functor<function_type, double( const double& )>::value == std::true_type , "function is of wrong type");
    if constexpr( is_functor<function_type, double( const double& )>::value )
    { 
      F.function = [](double x, void* p) -> double { return (*static_cast<function_type*>(p))(x); };
    }
    else if constexpr( is_functor<function_type, double( const std::array<double, 1>& )>::value )
    {
      F.function = [](double x, void* p) -> double { return (*static_cast<function_type*>(p))( std::array<double,1>{x} ); };
    }
    else static_assert( true, "function matches no signature!");

    F.params   = const_cast<function_type*>(&fcn);
    gsl_integration_qags (&F, min, max, 0, 1e-5, 1000, w, &result, &error);
    if( ws == nullptr ) gsl_integration_workspace_free (w);
    return result;
  }
  template <typename function_type> double integrate_1d_inf( const function_type& fcn, gsl_integration_workspace* ws = nullptr)
  {
    gsl_integration_workspace* w = (ws == nullptr) ? gsl_integration_workspace_alloc (1000) : ws;
    double result = 0;
    double error  = 0;
    gsl_function F;
    F.function = [](double x, void* p) -> double { return (*static_cast<function_type*>(p))(x); };
    F.params   = const_cast<function_type*>(&fcn);
    gsl_integration_qagi(&F, 0, 1e-5, 1000, w, &result, &error);
    if( ws == nullptr ) gsl_integration_workspace_free (w);
    std::cout << result << " +/- " << error << std::endl; 
    return result;
  }

  template <typename function_type> double integrate_1d_cauchy( const function_type& fcn, const double& x0, const double& min, const double& max, gsl_integration_workspace* ws = nullptr)
  {
    gsl_integration_workspace* w = (ws == nullptr) ? gsl_integration_workspace_alloc (1000) : ws;
    double result = 0;
    double error  = 0;
    gsl_function F;
    F.function = [](double x, void* p) -> double { return (*static_cast<function_type*>(p))(x); };
    F.params   = const_cast<function_type*>(&fcn);
    gsl_integration_qawc(&F, min, max, x0, 0, 1e-5, 1000, w, &result, &error);
    if( ws == nullptr ) gsl_integration_workspace_free (w);
    return result;
  }

  template <unsigned dim> struct integral {
    double value = {0};
    double var   = {0};
    int index = {0};
    std::array<double, dim> a     = {0};
    std::array<double, dim> b     = {0};
    bool operator<( const integral<dim>& other) const { return var < other.var; }
  };  

  template<unsigned dim, typename fcn> std::tuple<double, double, unsigned> integrate_fp( const fcn& F, const std::array<double,dim>& ctr, const std::array<double, dim>& wth )
  {
    if constexpr( is_functor<fcn, double(const std::array<float_v, dim>&)>::value and utils::size<float_v>::value == 4 )
      return AmpGen::integrate_fp_avx2d<dim>(F, ctr, wth);
    else 
    {
      return integrate_fp_scalar<dim>(F, ctr, wth );
    }
  }

  template <unsigned dim, typename fcn> double integrate(const fcn& F, const std::array<double, dim> xmin, const std::array<double, dim>&  xmax)
  {
    // Based off of ROOT Math/IntegratorMultiDim
    // With improved speed and safety:
    // - No dynamic memory allocation
    // - Static dimension of integrals 
    // - Supports vectorised integrands  
    // References to actual method [again, from ROOT Math/IntegratorMultiDim]
    //   1.A.C. Genz and A.A. Malik, Remarks on algorithm 006:
    //     An adaptive algorithm for numerical integration over
    //     an N-dimensional rectangular region, J. Comput. Appl. Math. 6 (1980) 295-302.
    //   2.A. van Doren and L. de Ridder, An adaptive algorithm for numerical
    //     integration over an n-dimensional cube, J.Comput. Appl. Math. 2 (1976) 207-217.
    if constexpr( dim == 1 ){
      if constexpr( is_functor<fcn, double(const double&)>::value ) { integrate_1d(F, xmin[0], xmax[0] ); }
      else if constexpr( is_functor<fcn, double(const std::array<double, 1>& ) >::value ){
        return integrate_1d( [&F](const double& x){ return F( std::array<double,1>{x} ); }, xmin[0], xmax[0] );
      }
      else static_assert(true, "1D function doesn't have recognised signature");
    }
    else {

      double epsrel = 1e-9;  //specified relative accuracy
      double epsabs = 0.; //specified relative accuracy
      //output parameters
      double relerr = 0 ; //an estimation of the relative accuracy of the result

      double result = 0;
      double abserr = 0;
      auto status  = 3;

      unsigned int ifncls = 0; /// number of function calls 
      bool  ldv   = false;
      unsigned isbrgn = 1;
      unsigned isbrgs = 1;

      constexpr unsigned irlcls = get_power<2,dim>::value +2*dim*(dim+1)+1; // number of function evaluations per iteration
      constexpr unsigned minpts = get_power<2,dim>::value +2*dim*(dim+1)+1; // minimum number of function evaluations
      constexpr unsigned maxpts = 100000;                                   // maximum number of function evaluations

      std::array<integral<dim>, ( 1 + maxpts/irlcls)/2 > partial_integrals; 

      integral<dim>* current_integral = &(partial_integrals[0]); 

      for (unsigned j=0; j<dim; j++) {
        current_integral->a[j] = (xmax[j] + xmin[j])*0.5;
        current_integral->b[j] = (xmax[j] - xmin[j])*0.5;
      }

      unsigned int idvax0=0; 
      integral<dim> tmp; 

      do {
        auto [rgnval, rgnerr, idvaxn] = integrate_fp<dim>( F, current_integral->a, current_integral->b );
        result += rgnval;
        abserr += rgnerr;
        ifncls += irlcls;

        double aresult = std::abs(result);
        tmp.var   = rgnerr;
        tmp.value = rgnval;
        tmp.index = idvaxn;
        tmp.a     = current_integral->a;
        tmp.b     = current_integral->b;

        if (ldv) {
          unsigned isbtmp =0;  
          while( true )
          {
            isbtmp = 2*isbrgn; 
            if (isbtmp > isbrgs) break;
            if (isbtmp < isbrgs && partial_integrals[isbtmp].var < partial_integrals[isbtmp+1].var ) isbtmp++;
            if (rgnerr >= partial_integrals[isbtmp].var ) break;
            partial_integrals[isbrgn] = partial_integrals[isbtmp] ;
            isbrgn = isbtmp;
          }
        }
        else {
          unsigned isbtmp =0;  
          do {
            isbtmp = isbrgn/2;
            if ( isbtmp >= 1 && rgnerr > partial_integrals[isbtmp].var ) {
              partial_integrals[isbrgn] = partial_integrals[isbtmp]; 
              isbrgn = isbtmp;
            }
          } while ( isbtmp >= 1 && rgnerr > partial_integrals[isbtmp].var );
        }

        partial_integrals[isbrgn] = tmp; 
        if (ldv) {//divison along chosen coordinate
          ldv = false; 
          current_integral = &(partial_integrals[isbrgn]);
          isbrgs += 1; 
          isbrgn  = isbrgs;
          partial_integrals[isbrgn].a     = current_integral->a;
          partial_integrals[isbrgn].b     = current_integral->b;
          partial_integrals[isbrgn].a[idvax0] += 2*current_integral->b[idvax0];
          current_integral = &(partial_integrals[isbrgn]);
          continue; 
        }
        //if no divisions to be made..
        relerr = std::abs(result) == 0 ? abserr : abserr/std::abs(result);
        if ((relerr < epsrel && aresult < epsabs) or 
            ( ( relerr < epsrel || abserr < epsabs ) && ifncls > minpts) ){ status = 0; break ; }
        if (( isbrgs >= partial_integrals.size()-1) or ( ifncls+2*irlcls > maxpts ) )  { status = 2; break; }
        //        std::cout << "#calls: " << ifncls << ", #cycles: " << isbrgs << " / " << ( 1 + maxpts/irlcls)/2 << std::endl; 
        ldv = true;
        isbrgn  = 1;
        current_integral = &( partial_integrals[isbrgn] );

        abserr -= current_integral->var;
        result -= current_integral->value;
        idvax0  = current_integral->index;

        current_integral -> b[idvax0] *= 0.5;
        current_integral -> a[idvax0] -= current_integral->b[idvax0];

      } while( status == 3 ); 

      return result;         //an approximate value of the integral

    }
  }
}
