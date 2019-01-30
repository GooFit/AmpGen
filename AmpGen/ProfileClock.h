#include <chrono>
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

namespace AmpGen{
  struct ProfileClock {
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start; 
    std::chrono::time_point<std::chrono::high_resolution_clock> t_end; 
    ProfileClock() : t_start(   std::chrono::high_resolution_clock::now()) {}
    void stop(){ t_end = std::chrono::high_resolution_clock::now() ; }
    operator double() const { return std::chrono::duration<double, std::milli>( t_end - t_start ).count() ; } ; 
  };

  template <int N, class FCN> 
    double  Profile( const FCN& fcn ){
      ProfileClock t; 
      for( size_t i = 0 ; i < N; ++i ) fcn();
      t.stop();
      INFO( typeof<FCN>() << " " << t/double(N) << "[ms] per iteration" );
      return t;
    }

}
