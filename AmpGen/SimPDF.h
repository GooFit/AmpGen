#ifndef AMPGEN_SIMPDF_H
#define AMPGEN_SIMPDF_H

#include <tuple>

namespace AmpGen
{
  class SimFit
  {
    std::vector<std::function<double( void )>> m_pdfs;

  public:
    SimFit() {}
    double getVal()
    {
      double LL = 0;
      for ( auto& pdf : m_pdfs ) LL += pdf();
      return LL;
    }

    template <class PDF>
    void add( PDF& pdf )
    {
      INFO( "Adding " << &pdf << " to sim fit" );
      m_pdfs.emplace_back( [&pdf]() -> double { return pdf.getVal(); } );
    }
  };
} // namespace AmpGen

#endif
