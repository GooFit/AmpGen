// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:17:55 GMT
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

template class AmpGen::NamedParameter<std::string>;
template class AmpGen::NamedParameter<double>;
template class AmpGen::NamedParameter<unsigned int>;
template class AmpGen::NamedParameter<int>;
template class AmpGen::NamedParameter<size_t>;
template class AmpGen::NamedParameter<bool>;
