#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/CombCorrLL.h"
#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/polyLASSO.h"
#include "AmpGen/ProfileClock.h"
#include <TMath.h>
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CombGamLL.h"
//#include <Math/IFunction.h>
#include <Math/Functor.h>
#include <TGraph.h>
#include <Minuit2/Minuit2Minimizer.h>
#include "AmpGen/PhaseCorrection.h"
#include <typeinfo>


#include <boost/algorithm/string/replace.hpp>
using namespace AmpGen;
using namespace std::complex_literals;
int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  size_t nFree =0 ;
  std::vector<std::string> freeParams;
  for (auto &p:MPS){
      if (p->isFree()){
          nFree++;
          freeParams.push_back(p->name());
      }
  }
  INFO("Have "<<nFree<<" free parameters");
  INFO("Fixing "<<freeParams[0]);
  MPS[freeParams[0]]->fix();
  nFree =0 ;
  for (auto &p:MPS){
      if (p->isFree()){
          nFree++;
          freeParams.push_back(p->name());
      }
  }
  INFO("Now have "<<nFree<<" free parameters");
  INFO(freeParams[0]<<" free? "<<MPS[freeParams[0]]->isFree());
 

  
  return 0;
}