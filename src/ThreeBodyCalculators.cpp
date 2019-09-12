#include <Math/AllIntegrationTypes.h>
#include <Math/IFunctionfwd.h>
#include <Math/ParamFunctor.h>
#include <Math/GSLIntegrator.h>
#include <Math/WrappedTF1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <memory.h>
#include <stddef.h>
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
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Projection.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/Units.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Event.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

template <class FCN> double dispersive( FCN& fcn , const double& s, double min , double max )
{
  TF1 fcn_tf1 = TF1( "fcn_tf1",fcn, min, max, 0 );
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0001);
  ig.SetFunction( ROOT::Math::WrappedTF1(fcn_tf1) );
  return ig.IntegralCauchy(min,max,s);
}

TGraph* ThreeBodyCalculator::runningMass(
    const double& mass,
    const double& min, 
    const double& max, 
    const size_t& nSteps, 
    const size_t& nSubtractions )
{
  double s0 = mass*mass;
  auto width0 = (*m_mps)[m_name+ "_width"]->mean() / getWidth(s0);
  auto fSubtracted = [&](const double* x, const double* p){return width0 * getWidth(*x)/pow(*x,nSubtractions);};
  double maxI = 1000;
  auto integral_term = [&](const double& s){
    return dispersive(fSubtracted,s,min,maxI) * mass * pow(s,nSubtractions) / M_PI;
  };  
  TGraph* dispersiveTerm = new TGraph();
  double g1 = integral_term( s0 - 0.001 );
  double g2 = integral_term( s0 + 0.001 );
  
  double A1 = (g2-g1)/(2*0.001);
  if( nSubtractions == 1 ) A1 = 0;
  double A0 = s0 - A1*s0 + integral_term(s0);
  if( nSubtractions == 0 ) A0 = 0;
  INFO("Calculated subtraction constants = [" << A0 << ", " << A1 << "]" );

  for ( size_t ind = 0; ind < nSteps; ++ind ) 
  {
    double s  = min + ( max - min ) * ind / double( nSteps );
    dispersiveTerm->SetPoint( ind, s,  A0 + A1 * s -integral_term(s)  );
    INFO( "Calculated dispersion relation as: " << A0 + A1 * s - integral_term(s) );
  }
  return dispersiveTerm;
}

TGraph* ThreeBodyCalculator::widthGraph( const size_t& steps, const double& min, const double& max )
{
  TGraph* g = new TGraph();
  double st = (max-min)/double(steps-1);
  for ( size_t c = 0; c < steps; ++c ) {
    double s = m_min + double( c ) * st;
    g->SetPoint( g->GetN(), s, getWidth(s) );
  }
  return g;
}

TGraph* ThreeBodyCalculator::fastRunningMass(
    const double& mass,
    const double& min, 
    const double& max, 
    const size_t& nSteps, 
    const size_t& nSubtractions )
{
  double maxI = 100 * GeV * GeV;
  double s0   = mass*mass;
  auto   g    = widthGraph(1000,min,maxI);
  INFO("Generated high-resolution width graph...");
  double gNorm = g->Eval(s0);
  double width = ParticlePropertiesList::get(m_name)->width();

  auto fSubtracted = [&](const double* x, const double* p){return g->Eval(*x)/(pow(*x,nSubtractions) ) ;};
  auto integral_term = [&](const double& s){
    return dispersive(fSubtracted,s,min,maxI) * width * mass * pow(s,nSubtractions) / ( gNorm * M_PI );
  };  
  TGraph* dispersiveTerm = new TGraph();
  double g1 = integral_term( s0 - 0.001 );
  double g2 = integral_term( s0 + 0.001 );
  
  double A1 = (g2-g1)/(2*0.001);
  if( nSubtractions == 1 ) A1 = 0;
  double A0 = s0 - A1*s0 + integral_term(s0);
  if( nSubtractions == 0 ) A0 = 0;
  for ( size_t ind = 0; ind < nSteps; ++ind ) 
  {
    double s  = min + ( max - min ) * ind / double( nSteps );
    dispersiveTerm->SetPoint( ind, s,  A0 + A1 * s -integral_term(s)  );
    INFO( "Calculated dispersion relation as: " << A0 + A1 * s - integral_term(s) );
  }
  return dispersiveTerm;
}


double ThreeBodyCalculator::PartialWidth::getWidth( const double& s )
{
  integrator.setMother( s );
  return integrator.integrate(totalWidth);
}

Expression ThreeBodyCalculator::PartialWidth::spinAverageMatrixElement(
    const std::vector<TransitionMatrix<complex_t>>& elements, DebugSymbols* msym )
{
  std::vector<Tensor> currents;
  for ( auto& element : elements ) {
    Particle particle(element.decayDescriptor(), type.finalStates() ); 
    auto perm = particle.identicalDaughterOrderings();
    for ( auto& p : perm ) {
      particle.setOrdering(p);
      particle.setLineshape( "FormFactor" );
      Expression prop = make_cse( element.coupling.to_expression() ) * make_cse( particle.propagator( msym ) );
      if ( msym != nullptr ) msym->emplace_back( element.decayTree.name() + "_g", element.coupling.to_expression() );
      if ( msym != nullptr ) msym->emplace_back( element.decayTree.name() + "_p", particle.propagator() );
      Tensor zt = particle.spinTensor(msym);
      zt.st() ;
      currents.push_back( zt * prop );
    }
  }
  Expression total;
  for ( auto& j_a : currents ) {
    for ( auto& j_b : currents ) total = total + dot( j_a, j_b.conjugate() );
  }
  ADD_DEBUG( total, msym );
  return pow(-1, ParticlePropertiesList::get(type.mother())->twoSpin() /2. ) * total;
}

ThreeBodyCalculator::ThreeBodyCalculator( const std::string& head, MinuitParameterSet& mps, const size_t& nKnots, const double& min, const double& max)
  : m_min(min),
    m_max(max),
    m_norm(1),
    m_nKnots(nKnots),
    m_name(head),
    m_mps(&mps)
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
    double s                   = m_min + double(c) * m_step;
    double I                   = getWidth(s);
    const std::string knotName = m_name + "::Spline::Gamma::" + std::to_string( c );
    if ( mps.find( knotName ) != nullptr ) mps[knotName]->setCurrentFitVal( I );
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
    double G = getWidth(s);
    INFO("Calculating width for " << c << " " << G );
    g->SetPoint( g->GetN(), s, G );
  }
  return g;
}

ThreeBodyCalculator::PartialWidth::PartialWidth( const EventType& evt, MinuitParameterSet& mps )
  : fcs( evt, mps, "" )
  , integrator(1, evt.mass(0)*evt.mass(0), evt.mass(1)*evt.mass(1) , evt.mass(2)*evt.mass(2) )
  , type(evt)
{
  DebugSymbols msym;
  Expression matrixElementTotal = spinAverageMatrixElement( fcs.matrixElements(), &msym );
  std::string name              = "";
  auto evtFormat = evt.getEventFormat();
  for ( auto& p : fcs.matrixElements() ) {
    name += p.decayDescriptor();
    partialWidths.emplace_back( spinAverageMatrixElement( {p}, &msym ), p.decayDescriptor(), evtFormat, DebugSymbols(), &mps );
  }
  totalWidth = CompiledExpression< std::complex<real_t>, const real_t*, const real_t* > ( matrixElementTotal, "width", evtFormat, {} , &mps );
  CompilerWrapper(true).compile( totalWidth, "");
}

double ThreeBodyCalculator::getWidth( const double& m )
{
  double G = 0;
  for ( auto& w : m_widths ){
    double wp = w.getWidth( m ); 
    if( std::isnan( wp) ) continue; 
    G += wp;
  }
  return G ; 
}

void ThreeBodyCalculator::setNorm( const double& mNorm )
{
  m_norm = 1;
  m_norm = getWidth(mNorm);
}

void ThreeBodyCalculator::makePlots(const double& mass, const size_t& x, const size_t& y)
{
  auto& sq      = m_widths[0].integrator;
  auto& evtType = m_widths[0].type;
  auto& fcs     = m_widths[0].totalWidth;
  auto projection_operators = evtType.defaultProjections( 500 );
  int points = NamedParameter<int>( "nPoints", 50000000 );
  sq.setMother( evtType.motherMass() );
  prepare();
  auto fcn = [&](const double* evt) { return std::real(fcs(evt)); };
  sq.makePlot( fcn, Projection2D( projection_operators[x], projection_operators[y] ), "s01_vs_s02", points )->Write();
}

void ThreeBodyCalculator::debug( const double& m, const double& theta )
{
  for( auto& width : m_widths )
  {
    Event event(12);
    width.integrator.setEvent({m,theta},event);
    event.print();
    width.totalWidth.debug( event.address() );
  }
}

