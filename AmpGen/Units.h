#ifndef AMPGEN_UNITS_H
#define AMPGEN_UNITS_H 1
#include "AmpGen/enum.h"

namespace AmpGen {
  
  static const double TeV = 1000;
  static const double GeV = 1;
  static const double MeV = 0.001;
  static const double KeV = 0.001*0.001;
  static const double  eV = 0.001*0.001*0.001; 

  declare_enum( Units, TeV, GeV, MeV, KeV, eV )
  double to_double(const Units& unit );
}
#endif
