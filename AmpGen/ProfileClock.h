#ifndef AMPGEN_PROFILECLOCK_H
#define AMPGEN_PROFILECLOCK_H 1
#include <chrono>
#include <math.h>
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

namespace AmpGen{
  struct ProfileClock {
    std::chrono::time_point<std::chrono::high_resolution_clock>       t_start;
    std::chrono::time_point<std::chrono::high_resolution_clock>       t_end;
    double t_duration = {0}; 
    
    ProfileClock() : t_start(std::chrono::high_resolution_clock::now()) {}
    void stop() 
    { 
      t_end      = std::chrono::high_resolution_clock::now() ; 
      t_duration += std::chrono::duration<double, std::milli>( t_end - t_start ).count();
    }
    double count() const 
    {
      auto now = std::chrono::high_resolution_clock::now() ; 
      return std::chrono::duration<double, std::milli>(now - t_start ).count();
    }
    void start(){ t_start = std::chrono::high_resolution_clock::now() ; }
    operator double() const { return t_duration; } ; 
  };

  template <int N, class FCN> 
    double  Profile( const FCN& fcn, const std::string& name ="" ){
      ProfileClock t; 
      for( size_t i = 0 ; i < N; ++i ) fcn();
      t.stop();
      INFO( (name == "" ? type_string<FCN>() : name ) << " " << t/double(N) << "[ms] per iteration" );
      return t;
    }
  template <int N, class FCN> 
    double  ProfileWithStat( const FCN& fcn, const std::string& name ="" ){
      double t = 0;
      double t2 = 0;
      double tmin = 1e9;
      double tmax = 0;
      for( size_t i = 0 ; i < N; ++i ){
        ProfileClock pi;
        fcn();
        pi.stop();
        t  += pi;
        t2 += pi*pi;
        tmin = pi < tmin ? pi : tmin;
        tmax = pi > tmax ? pi : tmax;
      }
      t /= double(N);
      t2 = std::sqrt( t2 / double(N) - t*t);
      INFO( (name == "" ? type_string<FCN>() : name ) << " " << t << " ± " << t2 << "[ms] per iteration << [" << tmin << ", " << tmax << "]" );
      return t;
    }
  
  template <int N, class FCN> 
    double  Profile2( const FCN& fcn ){
      ProfileClock t;
      auto z = 0 ; 
      for( size_t i = 0 ; i < N; ++i ) z += fcn();
      t.stop();
      INFO( type_string<FCN>() << " " << t/double(N) << "[ms] per iteration; " << z );
      return t;
    }
}
#endif
