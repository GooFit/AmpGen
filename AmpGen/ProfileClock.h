#include <chrono>

namespace AmpGen{
  struct ProfileClock {
    std::chrono::_V2::system_clock::time_point t_start; 
    std::chrono::_V2::system_clock::time_point t_end; 
    ProfileClock() : t_start(   std::chrono::high_resolution_clock::now()) {}
    void stop(){ t_end = std::chrono::high_resolution_clock::now() ; }
    operator double() const { return std::chrono::duration<double, std::milli>( t_end - t_start ).count() ; } ; 
  };
}
