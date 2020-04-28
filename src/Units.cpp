
#include "AmpGen/Units.h"
#include "AmpGen/Utilities.h"

namespace AmpGen {
  complete_enum(Units, TeV, GeV, MeV, KeV, eV)
}

double AmpGen::to_double(const AmpGen::Units& unit)
{
  static const double value_table[5] = {TeV, GeV, MeV, KeV, eV};
  return value_table[unsigned(unit)];  
} 
