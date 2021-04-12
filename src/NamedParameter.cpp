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

std::string AmpGen::optionalHelpString(const std::string& header, const std::vector<std::pair<std::string, std::string>>& args )
{
  std::string rt=header +"\n";
  for( size_t i = 0 ; i < args.size(); ++i ){
    rt += "\033[3m " + args[i].first  + std::string( 25 - args[i].first.size(), ' ');
    rt += "\033[0m: " + args[i].second + (i==args.size()-1 ? "" : "\n" );
  }
  return rt;
}
