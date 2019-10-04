#include "AmpGen/DTEventList.h"
using namespace std;
using namespace AmpGen;
#ifndef DTYIELDCALCULATOR
#define DTYIELDCALCULATOR
class DTYieldCalculator {
  public:
    DTYieldCalculator(const double& productionCrossSection = 3260) : 
      productionCrossSection(productionCrossSection){}
    double operator()(const double& lumi, 
        const AmpGen::EventType& t_signal, 
        const AmpGen::EventType& t_tag, 
        const bool& print = false);
    double bf( const AmpGen::EventType& type ) const;
  private: 
    double productionCrossSection;
    std::map<std::string, double> getKeyed( const std::string& name );
    std::map<std::string, double> branchingRatios = {getKeyed("BranchingRatios")};
    std::map<std::string, double> efficiencies    = {getKeyed("Efficiencies")};
};

//std::map<std::string, double> getKeyed(const std::string& name){
  //   std::map<std::string, double> s(name, 0.0);
   //  return s;
//}
std::map<std::string, double> DTYieldCalculator::getKeyed(const std::string& name)
{
  std::vector<std::string> things = AmpGen::NamedParameter<std::string>(name).getVector();
  std::map< std::string , double > branchingRatios; 
  for( auto& thing : things ){
    auto tokens = AmpGen::split( thing, ' ' );
    AmpGen::Particle p(tokens[0]);
    branchingRatios[p.uniqueString()]        = stod(tokens[1]);
    branchingRatios[p.conj().uniqueString()] = stod(tokens[1]);
  }
  return branchingRatios;
}

double DTYieldCalculator::operator()(const double& lumi, 
    const AmpGen::EventType& t_signal, 
    const AmpGen::EventType& t_tag, 
    const bool& print){
  auto statisticalFactor = 2;
  if( t_signal == t_tag || t_signal == t_tag.conj(false,true) ) statisticalFactor = 1;
  auto signal       = AmpGen::Particle(t_signal.decayDescriptor()).uniqueString(); 
  auto tag          = AmpGen::Particle(t_tag.decayDescriptor()).uniqueString();
  auto signalBar    = AmpGen::replaceAll( signal, "D0", "Dbar0");
  auto tagBar       = AmpGen::replaceAll( tag,    "D0", "Dbar0");
  auto eff          = [this](const std::string& tag) -> double { auto it = efficiencies.find(tag); 
    if( it == efficiencies.end() ){
      WARNING("Efficiency for final state: " << tag << " not found");
      return 1;
    }
    return it->second;
  };
  double efficiency = eff(signal) * eff(tag);
  double br         = branchingRatios[signal] * branchingRatios[tagBar] + branchingRatios[signalBar] * branchingRatios[tag];
  double totalDDbar = lumi * productionCrossSection; 
  if( print ){
    INFO("Expected yield for final state: " << t_signal << " vs " << t_tag );
    INFO("Total DDbar = " << totalDDbar );
    INFO("Efficiency  = " << efficiency );
    INFO("BR          = " << branchingRatios[signal] <<" x " << branchingRatios[tagBar] 
        << " + " << branchingRatios[signalBar] << " x " <<  branchingRatios[tag] 
        << " = " << br );
    INFO("Total       = " << statisticalFactor * totalDDbar * efficiency * br );
  }
  return statisticalFactor * totalDDbar * efficiency * br; 
}

double DTYieldCalculator::bf( const AmpGen::EventType& type ) const { 
  auto  p = AmpGen::Particle( type.decayDescriptor() );
  auto it = branchingRatios.find(p.decayDescriptor());
  if( it != branchingRatios.end() ) return it->second; 
  else {
    ERROR("Tag: " << p.decayDescriptor() << " not found in branching ratios");
    return 0;
  }
}


#endif