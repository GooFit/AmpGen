#include "AmpGen/Generator.h"
#include <TRandom3.h>

using namespace AmpGen;

extern "C" void AmpGen::python__generate(const char* eventType, double* out, const unsigned int size)
{
  auto rndm = TRandom3(1); /// FIXME: SEED is hardcoded 
  EventType type( split( std::string(eventType),' ') ); 
  auto phsp = Generator<PhaseSpace>(type, &rndm );
  auto events = phsp.generate( size ); 
  for( size_t i = 0 ; i < events.size(); ++i ){
    for( size_t j = 0 ; j < events[i].size(); ++j)
      out[events[i].size() * i + j] = events[i][j]; 
  }
}

