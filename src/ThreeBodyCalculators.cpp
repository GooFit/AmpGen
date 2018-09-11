#include <Math/AllIntegrationTypes.h>
#include <Math/IFunctionfwd.h>
#include <Math/ParamFunctor.h>
#include <TH2.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/DalitzIntegrator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/FastCoherentSum.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"

using namespace AmpGen;

TGraph* ThreeBodyCalculator::runningMass( const double& min, const double& max, const unsigned int& nSteps,
                                          const double& norm )
{

  setNorm( norm );

  TF1 kk = TF1( "kK",
                [&]( double* x, double* p ) {
                  double s = *x;
                  return getWidth( sqrt( s ) * 1000 );
                },
                min, 10000, 0 );
  ROOT::Math::WrappedTF1* wf1 = new ROOT::Math::WrappedTF1( kk );
  ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE );
  ROOT::Math::IGenFunction& stupid = *( wf1->Clone() );
  ig.SetFunction( stupid );
  ig.SetRelTolerance( 0.001 );
  TGraph* dispersiveTerm = new TGraph();
  for ( unsigned int ind = 0; ind < nSteps; ++ind ) {
    double s  = min + ( max - min ) * ind / double( nSteps );
    double ms = ig.IntegralCauchy( min, 10000, s );
    INFO( "Im(D[s=" << s << " = " << s << "]) = " << ms );
    dispersiveTerm->SetPoint( ind, s, ms );
  }
  return dispersiveTerm;
}

Tensor ThreeBodyCalculator::PartialWidth::zipTensor( const Tensor& tensor )
{
  Tensor zipped_tensor( tensor.dims() );
  for ( unsigned int i = 0; i < zipped_tensor.size(); ++i ) zipped_tensor[i] = SubTree( tensor[i] );
  return zipped_tensor;
}

double ThreeBodyCalculator::PartialWidth::getWidth( const double& m )
{
  integrator.setMother( m );
  double G = integrator.integrate( totalWidth );
  return G;
}

Expression ThreeBodyCalculator::PartialWidth::spinAverageMatrixElement(
    const std::vector<TransitionMatrix<std::complex<real_t>>>& elements, DebugSymbols* msym )
{
  struct Current {
    Expression propagator;
    Tensor spin_current;
  };
  std::vector<Current> currents;

  for ( auto& element : elements ) {
    auto perm = element.decayTree->identicalDaughterOrderings();
    for ( auto& p : perm ) {
      element.decayTree->setOrdering( p );
      Current c;
      std::string current_lineshape = element.decayTree->lineshape();
      element.decayTree->setLineshape( "FormFactor" );

      c.propagator = SubTree( element.coupling.to_expression() ) * SubTree( element.decayTree->Lineshape( msym ) );
      if ( msym != nullptr ) msym->emplace_back( element.decayTree->name() + "_g", element.coupling.to_expression() );
      if ( msym != nullptr ) msym->emplace_back( element.decayTree->name() + "_p", element.decayTree->Lineshape() );
      c.spin_current = zipTensor( element.decayTree->SpinTensor( msym ) );
      currents.push_back( c );
      element.decayTree->setLineshape( current_lineshape );
    }
  }
  Expression total;
  for ( auto& j_a : currents ) {
    if ( j_a.spin_current.size() == 4 ) {
      ADD_DEBUG( j_a.spin_current[0], msym );
      ADD_DEBUG( j_a.spin_current[1], msym );
      ADD_DEBUG( j_a.spin_current[2], msym );
      ADD_DEBUG( j_a.spin_current[3], msym );
    }
    for ( auto& j_b : currents ) {
      ADD_DEBUG( dot(j_a.spin_current,j_b.spin_current), msym );
      ADD_DEBUG( j_a.propagator * fcn::conj( j_b.propagator ) , msym );
      auto term = j_a.propagator * fcn::conj( j_b.propagator ) * dot( j_a.spin_current, j_b.spin_current );
      total     = total + term;
    }
  }
  ADD_DEBUG( total, msym );
  return -total;
}

ThreeBodyCalculator::ThreeBodyCalculator( const std::string& head, MinuitParameterSet& mps, const size_t& nKnots, const double& min, const double& max)
{
  std::vector<EventType> finalStates;
  auto rules                = AmplitudeRules( mps );
  auto rulesForThisParticle = rules.rulesForDecay( head );
  for ( auto& pAmp : rulesForThisParticle ) {
    auto type = pAmp.eventType();
    if ( std::find( finalStates.begin(), finalStates.end(), type ) != finalStates.end() ) continue;
    bool stateToBeExpanded = false;
    for ( auto& fs : type.finalStates() ) {
      if ( rules.hasDecay( fs ) ) stateToBeExpanded = true;
    }
    if ( stateToBeExpanded ) continue;
    INFO( "Final state to sum = " << type );
    finalStates.push_back( type );
  }
  for ( auto& type : finalStates ) m_widths.emplace_back( type, mps );

  bool isReady = true; 
  for( auto& width : m_widths ) isReady &= width.totalWidth.isReady();
  
  m_name                   = head;
  if( nKnots != 999) setAxis( nKnots, min, max ); 
}

void ThreeBodyCalculator::setAxis( const size_t& nKnots, const double& min, const double& max )
{
  m_nKnots = nKnots; 
  m_min    = min;
  m_max    = max;
  m_step   = ( m_max - m_min ) / double(m_nKnots-1);
}

void ThreeBodyCalculator::updateRunningWidth( MinuitParameterSet& mps, const double& mNorm )
{
  prepare();
  setNorm( mNorm == 0 ? mps[m_name + "_mass"]->mean() : mNorm );
  
  for ( size_t c = 0; c < m_nKnots; ++c ) {
    double s                   = m_min + double( c ) * m_step;
    double I                   = getWidth( sqrt( s ) * 1000 );
    const std::string knotName = m_name + "::Spline::Gamma::" + std::to_string( c );
    if ( mps.map().find( knotName ) != mps.map().end() ) mps[knotName]->setCurrentFitVal( I );
    INFO( knotName << " = " << I );
  }
}

void ThreeBodyCalculator::prepare()
{
  for ( auto& w : m_widths ) w.totalWidth.prepare();
}

TGraph* ThreeBodyCalculator::widthGraph( const double& mNorm )
{
  setNorm( mNorm );
  TGraph* g = new TGraph();
  for ( size_t c = 0; c < m_nKnots; ++c ) {
    double s = m_min + double( c ) * m_step;
    g->SetPoint( g->GetN(), s, getWidth( sqrt( s ) * 1000 ) );
  }
  return g;
}

ThreeBodyCalculator::PartialWidth::PartialWidth( const EventType& evt, MinuitParameterSet& mps )
    :
    fcs( evt, mps, "" )
    , integrator( 1000, evt.mass( 0 ), evt.mass( 1 ), evt.mass( 2 ) )
    , type( evt )
{
  DebugSymbols msym;
  Expression matrixElementTotal = spinAverageMatrixElement( fcs.matrixElements(), &msym );
  std::string name              = "";
  auto evtFormat = evt.getEventFormat();

  for ( auto& p : fcs.matrixElements() ) {
    name += p.decayTree->uniqueString();
    partialWidths.emplace_back( spinAverageMatrixElement( {p}, &msym ), p.decayTree->uniqueString(), evtFormat, DebugSymbols(), &mps );
  }
  totalWidth = CompiledExpression< std::complex<real_t>, const real_t*, const real_t* > ( matrixElementTotal, name, evtFormat, msym, &mps );
}

double ThreeBodyCalculator::getWidth( const double& m )
{
  double G = 0;
  for ( auto& w : m_widths ){
    double wp = w.getWidth( m );
    
    if( std::isnan( wp) ) continue; 
    
    G += wp;
  }
  return G ; // / m_norm;
}

void ThreeBodyCalculator::setNorm( const double& mNorm )
{
  m_norm = 1;
  m_norm = getWidth( mNorm );
}

void ThreeBodyCalculator::makePlots(const double& mass)
{

  auto& sq      = m_widths[0].integrator;
  auto& evtType = m_widths[0].type;
  if( mass != -1 ) evtType.setMotherMass( mass ); 
  auto& fcs     = m_widths[0].totalWidth;

  auto projection_operators = evtType.defaultProjections( 500 );
  int points = NamedParameter<int>( "nPoints", 50000000 );

  sq.setRandom( new TRandom3() );
  sq.setMother( evtType.motherMass() );

  prepare();

  std::function<double( const double* evt )> fcn = [&]( const double* evt ) {
    return std::real( fcs( evt ) );
  };

  sq.makePlot( fcn, Projection2D( projection_operators[0], projection_operators[1] ), "s01_vs_s02", points )->Write();

//  Projection sqDp1( [&]( const Event& evt ) { return sq.sqDp1( evt ); }, "m", "m^{'}", 500, 0, 1 );
//  Projection sqDp2( [&]( const Event& evt ) { return sq.sqDp2( evt ); }, "theta", "#theta^{'}", 500, 0, 1 );

//  sq.makePlot( fcn, Projection2D( sqDp1, sqDp2 ), "m_vs_theta", points )->Write();
}

void ThreeBodyCalculator::debug( const double& m, const double& theta )
{
  for( auto& width : m_widths )
  {
    Event event(12);
    width.integrator.setEvent({m,theta},event);
    event.print();
    width.totalWidth.debug( event );
  }
}

