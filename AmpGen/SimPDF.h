#ifndef AMPGEN_SIMPDF_H
#define AMPGEN_SIMPDF_H

#include <tuple>

namespace AmpGen
{
  class SimFit
  {
    public:
      SimFit() = default;
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
    
    private:
      std::vector<std::function<double( void )>> m_pdfs;

  };
} // namespace AmpGen

#endif
