#include <cmath>
#include <memory>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"

using namespace AmpGen;

DEFINE_LINESHAPE( Isotensor )
{
  //// I=2 pipi scattering
  //// Implements https://journals.aps.org/prd/pdf/10.1103/PhysRevD.78.052001
  //// we fix inelasticity to be zero, which should hold up to ~ 1.5 GeV

  Expression mpi( 0.139570 );
  Expression a          = Parameter( "Isotensor::a", 55.21 ) * M_PI / 180;
  Expression b          = Parameter( "Isotensor::b", 0.853 );
  Expression c          = Parameter( "Isotensor::c", -0.959 );
  Expression d          = Parameter( "Isotensor::d", 0.314 );
  Expression polyTerm   = pol( s, {1, b, c, d} );
  Expression phaseShift = -2 * a * fcn::sqrt( s / 4 - mpi * mpi ) / polyTerm;
  ADD_DEBUG( phaseShift, dbexpressions );
  ADD_DEBUG( fcn::sqrt( s / 4 - mpi * mpi ), dbexpressions );
  ADD_DEBUG( polyTerm, dbexpressions );
  Expression cos_shift = Cos( phaseShift );
  Expression sin_shift = Sin( phaseShift );
  Expression J =  Constant(0,1);
  return ( cos_shift + J*sin_shift - 1 ) * (-0.5 * J );
}
