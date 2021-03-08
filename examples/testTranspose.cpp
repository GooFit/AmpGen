#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Generator.h"
#include "TRandom3.h"
using namespace AmpGen;

int main(int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );
    MinuitParameterSet MPS;
    MPS.loadFromStream();
    TRandom3 rndm;
    rndm.SetSeed( 0 );
    gRandom = &rndm;

    EventType type = EventType(std::vector<std::string>({"D0", "K0S0", "pi-", "pi+"}));
    EventList data = Generator<>(type, &rndm).generate(1000);
    EventList mc = Generator<>(type, &rndm).generate(100000);

    CoherentSum A = CoherentSum(type, MPS);
    
    A.setEvents(data);
    A.setMC(mc);
    A.prepare();



    Event evt = data[0];
    Event evtT = evt;
    evtT.swap(1,2);


    INFO("A("<<evt.s(0,1)<<", "<<evt.s(0,2)<<") = "<<A.getVal(evt));
    INFO("A("<<evtT.s(0,1)<<", "<<evtT.s(0,2)<<") = "<<A.getVal(evtT));
    return 0;
}