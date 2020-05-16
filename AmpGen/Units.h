#ifndef AMPGEN_UNITS_H
#define AMPGEN_UNITS_H 1
#include "AmpGen/enum.h"

namespace AmpGen {
  
  static const double TeV = 1000;
  static const double GeV = 1;
  static const double MeV = 0.001;
  static const double KeV = 0.001*0.001;
  static const double  eV = 0.001*0.001*0.001; 

  static const double  ms = 1000*1000; 
  static const double  us = 1000; 
  static const double  ns = 1; 
  static const double  ps = 0.001; 
  static const double  fs = 0.001*0.001;

  static const double mm = 1.0;
  static const double um = 0.001; 
  static const double nm = 0.001*0.001; 
  declare_enum( Units, TeV, GeV, MeV, KeV, eV, ms, us, ns, ps, fs )
  double to_double(const Units& unit );
}
#endif
