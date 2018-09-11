#ifndef AMPGEN_MULTIQUARKCONTENT_H
#define AMPGEN_MULTIQUARKCONTENT_H
#include <iostream>
#include <string>
// Added to handle quark contents of mixed particles.
// For most particles, it will have only a single
// entry. But for mixed particles such as K-short
// with a quark content of p(dS)+q(Ds),
// MultiQuarkContent will contain one entry
// for each flavour option (i.e. dS and Ds, in this case).
// When adding up quark contents, all possibilities
// are considered (so if you add together the quark contents
// of two mixed particles, you get a MultiQuarkContent
// vector that contains four options).
// Used to see if a decay is comparible with being
// a strong/QED decay (quark contents compatible)
// or not. BW_BW.C uses this, for example, to
// determine if it should consider Parity conservation
// when determining the lowest possible angular
// momentum of a decay.
#include <vector>

#include "AmpGen/QuarkContent.h"

namespace AmpGen
{
  class MultiQuarkContent : public std::vector<QuarkContent>
  {
  public:
    MultiQuarkContent();
    MultiQuarkContent( const MultiQuarkContent& other );
    void antiThis();
    bool initFromString( const std::string& str );

    bool compatible( const MultiQuarkContent& other ) const;
    void print( std::ostream& os = std::cout ) const;

    MultiQuarkContent& operator+=( const MultiQuarkContent& rhs );
    MultiQuarkContent& operator-=( const MultiQuarkContent& rhs );

    MultiQuarkContent operator+( const MultiQuarkContent& rhs ) const;
    MultiQuarkContent operator-( const MultiQuarkContent& rhs ) const;

    bool operator==( const MultiQuarkContent& rhs ) const;
  };
  std::ostream& operator<<( std::ostream& st, const MultiQuarkContent& qc );
} // namespace AmpGen

#endif
//
