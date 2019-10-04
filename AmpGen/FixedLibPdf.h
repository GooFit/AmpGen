

#include "AmpGen/EventType.h"
#include "AmpGen/NamedParameter.h"

#include "AmpGen/OptionsParser.h"


#include "AmpGen/MinuitParameterSet.h"


#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


#include "AmpGen/EventList.h"




#include "AmpGen/MsgService.h"

#include "AmpGen/Utilities.h"






#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif
//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;

//#include <boost/algorithm/string.hpp>

class FixedLibPdf 
{
  public: 
    FixedLibPdf() = default; 
    FixedLibPdf(const EventType& type, MinuitParameterSet& mps) : 
      FixedLibPdf(NamedParameter<std::string>(type.decayDescriptor()+"::lib").getVal()) 
    {
      INFO("Constructing: " << type << " flib = " <<type.decayDescriptor()+"::lib" );
    }
    FixedLibPdf(const std::string& lib)
    {
      void* handle = dlopen( lib.c_str(), RTLD_NOW );
      if ( handle == nullptr ) ERROR( dlerror() );
      INFO("Constructing from " << lib );
      amp = AmpGen::DynamicFCN<complex_t( const double*, const int& )>( handle, "AMP" );
    }
    void prepare(){};
    void setEvents( AmpGen::EventList& evts ){};
    double prob_unnormalised( const AmpGen::Event& evt ) const { return std::norm(getValNoCache(evt)); }
    complex_t getValNoCache( const AmpGen::Event& evt ) const { return amp(evt,+1); }
    size_t size() { return 0; }
    void reset( const bool& flag = false ){};
  private:
    AmpGen::DynamicFCN<complex_t( const double*, const int& )> amp;
};