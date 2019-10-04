#ifndef AMPGEN_PROFILECLOCK_H
#define AMPGEN_PROFILECLOCK_H 1
#include <chrono>
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
    void start(){ t_start = std::chrono::high_resolution_clock::now() ; }
    operator double() const { return t_duration; } ; 
  };

  template <int N, class FCN> 
    double  Profile( const FCN& fcn ){
      ProfileClock t; 
      for( size_t i = 0 ; i < N; ++i ) fcn();
      t.stop();
      INFO( typeof<FCN>() << " " << t/double(N) << "[ms] per iteration" );
      return t;
    }
  
  template <int N, class FCN> 
    double  Profile2( const FCN& fcn ){
      ProfileClock t;
      auto z = 0 ; 
      for( size_t i = 0 ; i < N; ++i ) z += fcn();
      t.stop();
      INFO( typeof<FCN>() << " " << t/double(N) << "[ms] per iteration; " << z );
      return t;
    }
}
#endif
