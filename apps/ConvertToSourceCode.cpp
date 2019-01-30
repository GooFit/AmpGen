#include <complex>
#include <fstream>
#include <string>
#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ThreadPool.h"
#include "TRandom3.h"

using namespace AmpGen;


template <class T> void generate_source(T& pdf, EventList& normEvents, const std::string& sourceFile, MinuitParameterSet& mps, const double& sf)
{
  bool normalise                      = NamedParameter<bool>("Normalise",true);
  std::string type                    = NamedParameter<std::string>( "Type", "CoherentSum" );

  double norm = 1; 
  if( normalise ){
    double pMax = 0 ;
    pdf.setEvents( normEvents );
    pdf.prepare();

    for ( auto& evt : normEvents ) {
      if( type == "PolarisedSum" ){ 
        double px, py, pz; 
        gRandom->Sphere(px,py,pz, gRandom->Uniform(0,1));
        mps["Px"]->setCurrentFitVal(px);
        mps["Py"]->setCurrentFitVal(py);
        mps["Pz"]->setCurrentFitVal(pz);
        pdf.transferParameters();
      }
      double n = pdf.prob_unnormalised( evt );
      //    INFO( "Pol-vector =" << px << " " << py << " " << pz << " " << n1 << " " << n );
      //    INFO("=========================");
      if ( n > pMax ) {
        pMax = n;
      }
    }
    norm = pMax * sf ; 
  INFO( "Making binary with " << pMax << " x safety factor = " << sf );
  }
  pdf.generateSourceCode( sourceFile, norm, true );
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );
  //  ThreadPool::nThreads = 1;
  std::vector<std::string> oEventType = NamedParameter<std::string>( "EventType" ).getVector();
  std::string sourceFile              = NamedParameter<std::string>( "Output" , "output.cpp" );
  std::string type                    = NamedParameter<std::string>( "Type", "CoherentSum" );
  std::string outputPS                = NamedParameter<std::string>( "OutputEvents", "" );
  unsigned int NormEvents             = NamedParameter<unsigned int>( "NormEvents", 1000000 );
  double safetyFactor                 = NamedParameter<double>( "SafefyFactor", 3 );

  EventType eventType( oEventType );

  AmpGen::MinuitParameterSet MPS; //
  MPS.loadFromStream();
  /// This is just to calculate the overall normalisation of the PDF
  Generator<PhaseSpace> phsp( eventType );
  TRandom3 rnd;

  gRandom = &rnd;
  phsp.setRandom( &rnd );

  EventList phspEvents( oEventType );
  phsp.fillEventListPhaseSpace( phspEvents, NormEvents );

  if( type == "CoherentSum" ){
    CoherentSum sig( eventType, MPS, "" );
    generate_source( sig, phspEvents, sourceFile, MPS, safetyFactor );
  }
  if( type == "PolarisedSum" ){
    PolarisedSum sig( eventType, MPS );
    generate_source( sig, phspEvents, sourceFile, MPS, safetyFactor );
  }
  if ( outputPS != "" ) {
    std::ofstream ofile( outputPS );
    ofile << "0x,0y,0z,0t,1x,1y,1z,1t,2x,2y,2z,2t,3x,3y,3z,3t\n";
    for ( auto& event : phspEvents ) {
      for ( size_t i = 0; i < event.size(); i++ ) ofile << ( i == 0 ? "" : "," ) << event[i];
      ofile << "\n";
    }
  }
}
