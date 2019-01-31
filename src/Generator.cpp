#include "AmpGen/Generator.h"
#include <TRandom3.h>

using namespace AmpGen;

extern "C" void AmpGen::PyGenerate(const char* eventType, double* out, const unsigned int size)
{
  EventType type( split( std::string(eventType),' ') ); 
  INFO( type << " generating: " );
  auto phsp = Generator<PhaseSpace>(type, new TRandom3() );
  auto events = phsp.generate( size,0 ); 
  for( size_t i = 0 ; i < events.size(); ++i ){
    for( size_t j = 0 ; j < events[i].size(); ++j)
      out[events[i].size() * i + j] = events[i][j]; 
  }
}
