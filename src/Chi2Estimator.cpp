#include "AmpGen/Chi2Estimator.h"

#include <memory>
#include <ostream>

#include "AmpGen/ArgumentPack.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Event.h"

using namespace AmpGen;

struct Moment {
  double x;
  double xx;
  double N;
  std::vector<double> values;
  Moment() : x( 0 ), xx( 0 ), N( 0 ) {}
  void add( const double& value )
  {
    x += value;
    xx += value * value;
    N++;
    values.push_back( value );
  }
  void rescale( const double& val )
  {
    x *= val;
    xx *= ( val * val );
  }
  double val() { return x; }
  double var() { return N == 0 ? 0 : xx; }
};


void Chi2Estimator::doChi2( const EventList& dataEvents, const EventList& mcEvents,
    const std::function<double( const Event& )>& fcn )
{
  std::vector<Moment> data( m_binning.size() );
  std::vector<Moment> mc( m_binning.size() );

  INFO( "Splitting: " << dataEvents.size() << " data " << mcEvents.size() << " amongst " << m_binning.size()
      << " bins" );
  unsigned int j           = 0;
  double total_data_weight = 0;
  double total_int_weight  = 0;
  for ( auto& d : dataEvents ) {
    if ( j % 100000 == 0 ) INFO( "Binned " << j << " events" );
    double w = d.weight();
    data[m_binning.getBinNumber( d )].add( d.weight() );
    total_data_weight += w;
    j++;
  }
  j = 0;
  for ( auto& evt : mcEvents ) {
    if ( j % 100000 == 0 ) INFO( "Binned " << j << " events" );
    double w = fcn( evt ) * evt.weight() / evt.genPdf();
    mc[m_binning.getBinNumber( evt )].add( w );
    total_int_weight += w;
    j++;
  }
  double chi2 = 0;

  for ( unsigned int i = 0; i < m_binning.size(); ++i ) {
    mc[i].rescale( total_data_weight / total_int_weight );
    double delta = data[i].val() - mc[i].val();
    double tChi2 = delta * delta / ( data[i].val() + mc[i].var() );
    // if( tChi2 > 100 ) continue ;
    chi2 += tChi2;
  }
  m_chi2  = chi2;
  m_nBins = m_binning.size();
}

Chi2Estimator::Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
    const std::function<double( const Event& )>& fcn, const unsigned int& minEvents )
  : m_binning( dataEvents, MinEvents( minEvents ), Dim( dataEvents.eventType().dof() ) )
{
  doChi2( dataEvents, mcEvents, fcn );
}

Chi2Estimator::Chi2Estimator( const EventList& dataEvents, const EventList& mcEvents,
    const std::function<double( const Event& )>& fcn, const std::string& filename )
  : m_binning( File( filename ) )
{
  doChi2( dataEvents, mcEvents, fcn );
}
