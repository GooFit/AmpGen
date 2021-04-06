
#include "AmpGen/Units.h"
#include "AmpGen/Utilities.h"

namespace AmpGen {
  complete_enum(Units, TeV, GeV, MeV, KeV, eV, ms, us, ns, ps, fs)
}

double AmpGen::to_double(const AmpGen::Units& unit)
{
  static const double value_table[10] = {TeV, GeV, MeV, KeV, eV, ms, us, ns, ps, fs};
  return value_table[unsigned(unit)];  
} 
