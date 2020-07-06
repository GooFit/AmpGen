#include "AmpGen/gCoherentSum.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ratio>
#include <thread>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/ThreadPool.h"
#include "AmpGen/ProfileClock.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace AmpGen;


gCoherentSum::gCoherentSum(const EventType& type, const MinuitParameterSet& mps):
    m_type(type),
    m_typeConj(type.conj(true))
{
    m_A = CoherentSum(m_type, mps);
    m_C = CoherentSum(m_typeConj, mps);
    m_x = mps["gCoherentSum::x"];
    m_y = mps["gCoherentSum::y"];



    INFO("Making gCoherentSum");
    prepare();

}
