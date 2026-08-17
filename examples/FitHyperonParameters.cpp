#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/Property.h"

using namespace AmpGen;

int main(int argc, char *argv[]) {
  OptionsParser::setArgs(argc, argv);
  MinuitParameterSet mps;
  mps.add("gRe", AmpGen::Flag::Free, 0, 0.1, -10, 10);
  mps.add("gIm", AmpGen::Flag::Free, 0, 0.1, -10, 10);
  const auto alpha_measured = Property<double>(nullptr, "alpha", 0.75);
  const auto phi_measured = Property<double>(nullptr, "phi", -6.5);
  const auto alpha_error = Property<double>(nullptr, "alpha_error", 0.012);
  const auto phi_error = Property<double>(nullptr, "phi_error", 3.5);
  std::function<double()> f = [&mps, &alpha_measured, &alpha_error, &phi_measured, &phi_error]() -> double {
    auto alpha_pred = 2. * mps[0]->mean() / (1. + mps[0]->mean() * mps[0]->mean() + mps[1]->mean() * mps[1]->mean());
    auto phi_pred = 180. * atan(2. * mps[1]->mean() / (1. - (mps[0]->mean() * mps[0]->mean() + mps[1]->mean() * mps[1]->mean()))) / M_PI;
    return pow((alpha_measured - alpha_pred) / alpha_error, 2.) + pow((phi_measured - phi_pred) / phi_error, 2.);
  };
  Minimiser(f, &mps).doFit();
}
