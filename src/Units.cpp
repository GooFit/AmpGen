
#include "AmpGen/Units.h"
#include "AmpGen/Utilities.h"

namespace AmpGen {
  complete_enum(Units, TeV, GeV, MeV, KeV)
}

double AmpGen::to_double(const AmpGen::Units& unit)
{
  static constexpr std::array<double,4> value_table = {TeV, GeV, MeV, KeV};
  return value_table[unsigned(unit)];  
} 
