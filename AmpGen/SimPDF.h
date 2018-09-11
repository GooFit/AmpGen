#ifndef AMPGEN_SIMPDF_H
#define AMPGEN_SIMPDF_H

#include <tuple>

namespace AmpGen
{
  /*
     template<std::size_t I = 0, typename FuncT, typename... Tp>
     inline typename std::enable_if<I == sizeof...(Tp), void>::type
     for_each(std::tuple<Tp...> &, FuncT) // Unused arguments are given no names.
     { }

     template<std::size_t I = 0, typename FuncT, typename... Tp>
     inline typename std::enable_if<I < sizeof...(Tp), void>::type
     for_each(std::tuple<Tp...>& t, FuncT f)
     {
     f(std::get<I>(t));
     for_each<I + 1, FuncT, Tp...>(t, f);
     }
     template < class ...TYPES>
     class SimLL {
     SimLL( const TYPES & ...  _pdfs ) : m_pdfs( std::tuple<TYPES...>(_pdfs...) ) {}
     std::tuple<TYPES...> m_pdfs;
     double getVal(){
     double LL=0;
     for_each( m_pdfs, [&LL](auto& f){ LL+=f.getVal() ; } );
     }
     }
     */
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
