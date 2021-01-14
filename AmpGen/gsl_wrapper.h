
#include <gsl/gsl_integration.h>

namespace AmpGen {
  template <typename function_type> double integrate_1d( const function_type& fcn, const double& min, const double& max, gsl_integration_workspace* ws = nullptr)
  {
    gsl_integration_workspace* w = (ws == nullptr) ? gsl_integration_workspace_alloc (1000) : ws;
    double result = 0;
    double error  = 0;
    gsl_function F;
    F.function = [](double x, void* p) -> double { return (*static_cast<function_type*>(p))(x); };
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


  template <typename function_type> double integrate_2d( const function_type& fcn, const double& x1, const double& x2, const double& y1, const double& y2)
  {
    gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000); 
    auto outer = [fcn, &x1, &x2, &w1](const double& y)
    {
      auto inner = [fcn, &y](const double& x){ return fcn(x,y); };
      return integrate_1d(inner, x1,x2, w1); 
    };
    auto r = integrate_1d( outer, y1, y2);
    gsl_integration_workspace_free (w1);
    return r;
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
}
