#include "AmpGen/FitResult.h" 

using namespace AmpGen; 

int main( int argc, char** argv )
{
  if( argc < 2 ) 
  {
    FATAL("Requires parameters fit result, [optional: input file]");
  }
  FitResult(argv[1]).writeOptions(std::cout, argc == 3 ? argv[2] : "" ); 
}
