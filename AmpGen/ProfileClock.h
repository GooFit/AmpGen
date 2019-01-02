#include <chrono>

namespace AmpGen{
  struct ProfileClock {
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start; 
    std::chrono::time_point<std::chrono::high_resolution_clock> t_end; 
    ProfileClock() : t_start(   std::chrono::high_resolution_clock::now()) {}
    void stop(){ t_end = std::chrono::high_resolution_clock::now() ; }
    operator double() const { return std::chrono::duration<double, std::milli>( t_end - t_start ).count() ; } ; 
  };
}
