#include <ctype.h>
#include <stdlib.h>
#include <memory.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <complex>

#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Vertex.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/QuarkContent.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Tensor.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Units.h"
#include "AmpGen/Wigner.h"
#include "AmpGen/Types.h"

using namespace AmpGen;
using AmpGen::Lineshape::LineshapeFactory;

Particle::Particle( const std::string& name, const Particle& p1, const Particle& p2 ) 
{
  m_props = ParticlePropertiesList::get( name );
  m_name  = name;
  m_daughters.push_back( std::make_shared<Particle>( p1 ) );
  m_daughters.push_back( std::make_shared<Particle>( p2 ) );
  pdgLookup();
  sortDaughters();
  for ( auto& d : m_daughters ) d->setTop( false );
  m_uniqueString = makeUniqueString();
}

Particle::Particle( const int& pdgID, const Particle& p1, const Particle& p2 )
{
  m_props = ParticlePropertiesList::get( pdgID );
  m_daughters.push_back( std::make_shared<Particle>( p1 ) );
  m_daughters.push_back( std::make_shared<Particle>( p2 ) );
  if ( m_props != nullptr ) m_name = m_props->name();
  pdgLookup();
  sortDaughters();
  for ( auto& d : m_daughters ) d->setTop( false );
  m_uniqueString = makeUniqueString();
}
Particle::Particle( const std::string& decayString, const std::vector<std::string>& finalStates,
    const bool& orderDaughters )
  : Particle()
{
  DEBUG( "Building particle from decay string = " << decayString );
  auto items = getItems( decayString );
  if ( items.size() == 0 ) {
    ERROR( "Failed to parse string" << decayString );
    m_isStateGood = false;
    return;
  }
  m_modifiers = getItems( items[0], {"[", "]"}, ";" );
  m_name      = m_modifiers[0];
  m_modifiers.erase( m_modifiers.begin() );
  for ( auto& modifier : m_modifiers ) parseModifier( modifier );

  for ( unsigned int it = 1; it < items.size(); ++it ) m_daughters.push_back( std::make_shared<Particle>( items[it] ) );

  m_props = ParticlePropertiesList::get( m_name );
  pdgLookup();
  auto fs = getFinalStateParticles( false );
  std::vector<bool> hasUsedFinalState( finalStates.size(), 0 );
  for ( auto d : fs ) {
    for ( unsigned int i = 0; i < finalStates.size(); ++i ) {
      if ( hasUsedFinalState[i] || d->name() != finalStates[i] ) continue;
      DEBUG( "HEAD = " << name() << " " << d << " name = " << d->name() << " setting index = " << i );
      d->setIndex( i, true );
      hasUsedFinalState[i] = true;
      break;
    }
  }
  if ( orderDaughters ) sortDaughters();
  for ( auto& d : m_daughters ) m_isStateGood &= d->isStateGood();

  if ( !isStateGood() ) {
    ERROR( "Amplitude " << decayString << " not configured correctly" );
  }
  if ( finalStates.size() == fs.size() ) {
    std::string finalStateString = "";
    for ( auto& fs : finalStates ) finalStateString += fs + " ";
    for ( auto used : hasUsedFinalState ) {
      m_isStateGood &= used;
    }
    if ( !isStateGood() ) {
      WARNING( "Amplitude " << decayString << " does not match requested event type " << finalStateString );
    }
  }
  m_uniqueString = makeUniqueString();
}

Particle::Particle( const std::string& name, const unsigned int& index ) : Particle()
{
  m_props = ParticlePropertiesList::get( name );
  m_index = index;
  pdgLookup();
  for ( auto& d : m_daughters ) d->setTop( false );
  m_uniqueString = makeUniqueString();
}

void Particle::parseModifier( const std::string& mod )
{
  if( mod.size() == 1 )
  {
    DEBUG( "Modifier = " << mod );
    if ( mod == "S" ) m_orbital = 0;
    else if ( mod == "P" ) m_orbital = 1;
    else if ( mod == "D" ) m_orbital = 2;
    else if ( mod == "F" ) m_orbital = 3;
    else if ( mod == "G" ) m_orbital = 4;
    else if ( mod == "H" ) m_orbital = 5;
    else if ( mod == "I" ) m_orbital = 6;
    bool status=true;
    auto tmp = lexical_cast<unsigned int>(mod,status);
    if( status ) m_spinConfigurationNumber = tmp;
  }
  else if ( mod.size() == 2 )
  {
    parseModifier( mod.substr(0,1) );
    parseModifier( mod.substr(1,1) );
  }
  else if ( LineshapeFactory::isLineshape( mod ) )
    m_lineshape = mod;
}

double Particle::spin() const { return double( m_props->twoSpin() / 2. ) ; }
double Particle::S() const { return m_spinConfigurationNumber ; }

std::string Particle::vertexName() const
{
  if ( m_daughters.size() != 2 ){
    WARNING("Vertices only well-defined for quasi two-body processes, check logic");
    return "ERROR";
  }
  std::string orbitalString; 
  if( m_orbital == 0 ) orbitalString = "S";
  if( m_orbital == 1 ) orbitalString = "P";
  if( m_orbital == 2 ) orbitalString = "D";
  if( m_orbital == 3 ) orbitalString = "F";
  if( m_spinConfigurationNumber != 0 ) orbitalString += std::to_string( m_spinConfigurationNumber );

  auto vx=  m_props->spinName() + "_" + 
    m_daughters[0]->m_props->spinName() + 
    m_daughters[1]->m_props->spinName() + "_" + orbitalString;

  return vx;
}

void Particle::pdgLookup()
{
  if ( m_props == nullptr ) {
    m_isStateGood = false;
    return;
  }

  if ( m_lineshape == "BW" || m_usesDefaultLineshape ) {
    if ( m_name.find( "NonRes" ) != std::string::npos || m_props->width() < 1e-3 ) 
      m_lineshape = "FormFactor";
    m_usesDefaultLineshape = true;
  } else
    m_usesDefaultLineshape = false;
  m_parity                 = m_props->P();
  if ( !isdigit( m_props->J()[0] ) ) {
    ERROR( "Spin not recognised! : " << m_name << " J = " << m_props->J() );
  }
  std::string default_modifier = NamedParameter<std::string>("Particle::DefaultModifier","");
  if( default_modifier != "" && m_lineshape.find(".") == std::string::npos ){
    m_lineshape = m_lineshape + "." +default_modifier;
  }
  bool isStrong = quarks() == daughterQuarks();
  DEBUG( m_name << " is decaying via " << ( isStrong ? "strong" : "electroweak" ) << " interactions; P = " << props()->P() );

  if ( m_name.find( "NonRes" ) != std::string::npos ) isStrong = true;
  m_minL = m_daughters.size() == 2 ? orbitalRange( isStrong ).first : 0;
  if ( m_daughters.size() == 2 ) {
    DEBUG( "IsStrong ? " << isStrong << " orbital =  " << m_orbital << " minL   =  " << m_minL
        << " name   =  " << m_name << " d0     =  " << m_daughters[0]->name()
        << " d1     =  " << m_daughters[1]->name() );
  }
  if ( m_orbital == 0 ) m_orbital = m_minL; //// define in ground state unless specified ///

  for ( auto& d : m_daughters ) d->setTop( false );
}

Tensor Particle::P() const
{
  Tensor momentum( std::vector<double>( {0., 0., 0., 0.} ), std::vector<size_t>( {4} ) );

  if ( isStable() ) {
    if ( m_index != 999 ) {
      const std::string index = std::to_string( m_index );
      Tensor rt( std::vector<Expression>( {
            Parameter( index + "_Px", 0, false, 1 ), 
            Parameter( index + "_Py", 0, false, 1 ),
            Parameter( index + "_Pz", 0, false, 1 ), 0 } ) , {4} );
      rt[3] = make_cse( fcn::sqrt( mass()*mass() + rt[0]*rt[0] + rt[1]*rt[1] + rt[2]*rt[2] ) );
         //   Parameter( index + "_E" , 0, false, 1 )} ),
      return rt;
    } else
      ERROR( "Stable particle " << m_index << "is unindexed!" );
  } else {
    for ( auto& d : m_daughters ) momentum = momentum + d->P();
  }
  return momentum;
}

Tensor Particle::Q() const
{
  if ( m_daughters.size() != 2 ) {
    ERROR( " Q is only well defined for 2 body decay - check particle logic" );
    return Tensor();
  }
  return m_daughters[0]->P() - m_daughters[1]->P();
}

std::shared_ptr<Particle> Particle::daughter( const size_t& index ) { return ( m_daughters[index] ); }
std::shared_ptr<Particle> Particle::daughter( const size_t& index ) const { return m_daughters[index]; }

std::string Particle::orbitalString() const
{
  std::string orbital_part; 
  if ( m_orbital == 0 ) orbital_part = "S";
  else if ( m_orbital == 1 ) orbital_part = "P";
  else if ( m_orbital == 2 ) orbital_part = "D";
  else if ( m_orbital == 3 ) orbital_part = "F";
  else orbital_part = "X";

  if( m_spinConfigurationNumber == 0 ) return orbital_part; 
  else return orbital_part + std::to_string( m_spinConfigurationNumber );
}

bool Particle::hasModifier( const std::string& modifier ) const
{
  return std::find( m_modifiers.begin(), m_modifiers.end(), modifier ) != m_modifiers.end();
}

std::string Particle::modifierString() const
{
  std::string modString = "[";
  auto append = [&modString](const std::string& toAppend){ 
    modString += ( modString == "[" ? toAppend : ";" +toAppend );
  };
  if ( m_orbital != m_minL || m_spinConfigurationNumber != 0 ) append( orbitalString() );
  if ( !m_usesDefaultLineshape ) append( m_lineshape );
  std::vector<std::string> otherModifiers = {"BgSpin0", "Inco", "NoSym"};
  for ( auto& mod : otherModifiers ) if ( hasModifier( mod ) ) append(mod) ;
  for ( auto& mod : m_modifiers ) if( mod.find("=") != std::string::npos ) append( mod );
  modString += "]";
  return modString == "[]" ? "" : modString;
}

std::string Particle::makeUniqueString()
{
  std::string modifier = modifierString();
  if ( m_daughters.size() != 0 ) {
    std::string val = m_name + modifier + "{";
    for ( unsigned int i = 0; i < m_daughters.size(); ++i ) {
      m_daughters[i]->makeUniqueString();
      val += m_daughters[i]->uniqueString() + ( i == m_daughters.size() - 1 ? "}" : "," );
    }
    m_uniqueString = val;
  } else {
    m_uniqueString = m_name + modifier ;
  }
  return m_uniqueString;
}

std::vector<std::shared_ptr<Particle>> Particle::getFinalStateParticles( const bool& sort ) const
{
  std::vector<std::shared_ptr<Particle>> ffs;

  for ( auto daughter : m_daughters ) {
    if ( daughter->isStable() ) {
      ffs.push_back( daughter );
    } else {
      auto daughter_states = daughter->getFinalStateParticles( sort );
      DEBUG( "Pushing back daughters of: " << name() );
      for ( auto& s : daughter_states ) {
        DEBUG( "Pushing back: " << s->name() << " " << s );
        ffs.push_back( s );
      }
    }
  }
  if ( sort ) {
    std::sort( ffs.begin(), ffs.end(), []( auto& p1, auto& p2 ) { return p1->originalIndex() < p2->originalIndex(); } );
  }
  return ffs;
}

Particle Particle::quasiStableTree() const
{
  Particle returnVal( m_name );
  for ( auto& d : m_daughters ) {
    auto qTree = d->quasiStableTree();
    if ( d->isQuasiStable() ) {
      returnVal.addDaughter( std::make_shared<Particle>( qTree ) );
    } else {
      for ( auto& dp : qTree.daughters() ) returnVal.addDaughter( dp );
    }
  }
  returnVal.sortDaughters();
  returnVal.makeUniqueString();
  return returnVal;
}

Expression Particle::propagator( DebugSymbols* db ) const
{
  if ( db != nullptr && !isStable() ) db->emplace_back( uniqueString() +" lineshape", Parameter( "NULL", 0, true ) );
  if ( m_daughters.size() == 0 ) return 1;

  Expression total( 1. );
  DEBUG( "Getting lineshape " << m_lineshape << " for " << m_name );
  //Expression s = m_isHead ? dot(P(),P()) : massSq();
  Expression s= massSq();
  if ( m_daughters.size() == 2 ) {
    total = total * make_cse( LineshapeFactory::getLineshape( m_lineshape, s, daughter( 0 )->massSq(), daughter( 1 )->massSq(), m_name, m_orbital, db ) );
  }
  else if ( !m_isHead && m_daughters.size() == 3 ) {
    DEBUG( "Three-body propagator defaults to fixed-width Breit-Wigner" );
    std::string shape = m_lineshape == "BW" ? "SBW" : m_lineshape;
    auto propagator   = ( LineshapeFactory::getGenericShape(
          shape, {daughter( 0 )->P(), daughter( 1 )->P(), daughter( 2 )->P()}, m_name, m_orbital, db ) );

    total = total * propagator ;  
  }  
  for ( auto& d : m_daughters ) total = total * make_cse( d->propagator( db ) );
  if ( db != nullptr ) db->emplace_back( "A(" + uniqueString() + ")", total );
  return total;
}

std::vector<std::vector<size_t>> Particle::identicalDaughterOrderings() const
{
  std::vector<std::vector<size_t>> orderings;
  auto finalStateParticles = getFinalStateParticles();
  std::vector<std::string> permutation0_names;
  std::vector<size_t> indices( finalStateParticles.size() );
  std::iota( indices.begin(), indices.end(), 0 );
  for ( auto& particle : finalStateParticles ) permutation0_names.push_back( particle->name() );
  do {
    bool isSameAsP0 = true;
    for ( size_t i = 0; i < indices.size(); ++i )
      if ( permutation0_names[i] != finalStateParticles[indices[i]]->name() ) isSameAsP0 = false;
    if ( isSameAsP0 ) orderings.push_back( indices );
  } while ( std::next_permutation( indices.begin(), indices.end() ) );

  return orderings;
}

void Particle::setOrdering( const std::vector<size_t>& ordering )
{
  auto finalStateParticles = getFinalStateParticles();
  for ( unsigned int i = 0; i < ordering.size(); ++i ) {
    finalStateParticles[i]->setIndex( ordering[i] );
  }
}
void Particle::addDaughter( const std::shared_ptr<Particle>& particle ) { m_daughters.push_back( particle ); }

Tensor Particle::transitionMatrix( DebugSymbols* db  )
{
  auto particles = getFinalStateParticles();
  return spinTensor();
}

Expression Particle::getExpression( DebugSymbols* db, const unsigned int& index )
{
  if( db != nullptr && !isStable() ) 
    db->emplace_back( uniqueString() , Parameter( "NULL", 0, true ) );

  std::string spinFormalism = NamedParameter<std::string>("Particle::SpinFormalism","Covariant");
  auto localFormalism = attribute("SpinFormalism"); 
  if( localFormalism.first ) spinFormalism = localFormalism.second;
  Expression total( 0 );
  auto finalStateParticles  = getFinalStateParticles();
  auto orderings            = identicalDaughterOrderings();
  bool doBoseSymmetrisation = !hasModifier( "NoSym" );
  bool sumAmplitudes        = !hasModifier( "Inco" );
  bool includeSpin          = !hasModifier( "BgSpin0" );
  for ( auto& ordering : orderings ) {
    setOrdering( ordering );
    Expression spinFactor = 1; 
    if( includeSpin && spinFormalism == "Covariant" ){
      Tensor st = spinTensor(db);
      if( m_props->twoSpin() == 0 ) spinFactor = st[0];
      if( m_props->twoSpin() == 1 ){
        Tensor is = Bar( externalSpinTensor(m_polState) );
        ADD_DEBUG_TENSOR( externalSpinTensor(m_polState), db );
        ADD_DEBUG_TENSOR( st, db );
        Tensor::Index a;
        if( st.size() != 4 ){ 
          ERROR("Spin tensor is the wrong rank = " << st.dimString() ); 
          spinFactor = 1; 
        }
        else { spinFactor = is(a) * st(a) ; }
      }
    }
    if ( includeSpin && spinFormalism == "Wigner" ){
      spinFactor = helicityAmplitude( *this, TransformSequence(), double(polState())/2.0, db ); 
    }
    if( db != nullptr ){
      std::string finalStateString="";
      for( auto& fs : finalStateParticles ) 
        if( fs->spin() != 0 ){
          finalStateString += std::to_string(fs->polState()) + "_";
        }
      if( finalStateString != "" ) 
        finalStateString = finalStateString.substr(0, finalStateString.size()-1);
      db->emplace_back( "SF_"+std::to_string(polState()) +"_"+ finalStateString , spinFactor );
    }
    DEBUG( "Got spin matrix element -> calculating lineshape product" );
    if ( sumAmplitudes ) total = total + make_cse ( propagator( db ) ) * spinFactor;
    else {
      Expression ls     = propagator( db ) * spinFactor;
      Expression conjLs = fcn::conj( ls );
      Expression prod   = ls * conjLs;
      ADD_DEBUG( ls, db );
      ADD_DEBUG( conjLs, db );
      ADD_DEBUG( prod, db );
      total = total + prod;
    }
    if ( !doBoseSymmetrisation ) break;
  }
  ADD_DEBUG( total, db );
  double nPermutations = doBoseSymmetrisation ? orderings.size() : 1;
  if ( sumAmplitudes ) return total / fcn::sqrt( nPermutations );
  else {
    Expression sqrted = fcn::sqrt( total / nPermutations );
    ADD_DEBUG( sqrted, db );
    return sqrted;
  }
}

Tensor Particle::spinTensor( DebugSymbols* db ) const
{
  DEBUG( "Getting SpinTensor for : " << m_name << " " << m_daughters.size() );
  if ( m_daughters.size() == 0 ){
    auto S = externalSpinTensor(m_polState, db);
    if( S.size() != 1 ) ADD_DEBUG_TENSOR_NAMED( S, db, "S("+name()+")" );
    return S;
  }
  else if ( m_daughters.size() == 2 ) {
    Tensor value = VertexFactory::getSpinFactor( P(), Q(), 
        daughter(0)->spinTensor( db ),
        daughter(1)->spinTensor( db ), vertexName() , db );
    DEBUG( "Returning spin tensor" );
    return value;
  } else if ( m_daughters.size() == 3 ) {
    return VertexFactory::getSpinFactorNBody( {
        {daughter(0)->P(), daughter(0)->spinTensor()},
        {daughter(1)->P(), daughter(1)->spinTensor()},
        {daughter(2)->P(), daughter(2)->spinTensor()}}, spin() , db );
  } else if ( m_daughters.size() == 1 ) {
    DEBUG( "Forwarding through Quasi-particle" );
    return daughter( 0 )->spinTensor( db );
  }
  return Tensor( std::vector<double>( {1.} ), {1} );
}

Tensor Particle::externalSpinTensor(const int& polState, DebugSymbols* db ) const
{
  DEBUG("Getting final state spin tensor for: " << name() << " " << spin() );
  if ( spin() == 0 )
    return Tensor( std::vector<double>( {1.} ), std::vector<size_t>( {1} ) );

  std::string basis = NamedParameter< std::string >("Particle::SpinBasis","Dirac");
  Tensor p        = P();
  Expression pX   = p.get(0) ;
  Expression pY   = p.get(1) ;
  Expression pZ   = p.get(2) ;
  Expression pE   = p.get(3) ;
  Expression pP   = fcn::sqrt( pX*pX + pY*pY + pZ*pZ );
  Expression m    = fcn::sqrt( dot( p, p ) );
  complex_t I(0,1);
  Expression z    = pX + I*pY; 
  Expression zb   = pX - I*pY;
  auto id = props()->pdgID(); 
  if ( m_props->twoSpin() == 2 && basis == "Weyl" ) {
    Expression N = 1 / ( m * fcn::sqrt( 2 ) );
    Expression em = pE+m;
    if( polState ==  0 ) return    Tensor( { pX*pZ/(m*em), pY*pZ/(m*em), 1. + pZ*pZ/(m*em), pZ/m } );
    if( polState ==  1 ) return -N*Tensor( { m + pX*zb/em, pY*zb/em -m*I, pZ*zb, zb } ); 
    if( polState == -1 ) return  N*Tensor( { m + pX*z/em, pY*z/em + m*I, pZ*z, z } ); 
  }
  if( m_props->twoSpin() == 2 && basis == "Dirac" ){
    Expression invPT2 = Ternary( pX*pX + pY*pY > 1e-6, 1./(pX*pX+pY*pY), 0 );
    if( m_props->mass() == 0  ){
      Expression N = 1./(pP*sqrt(2)); 
      if( polState ==  1 ) return N * Tensor({ -pZ - pY*invPT2*(pP-pZ)*(pY-I*pX), 
          -I*pZ + pX*invPT2*(pP-pZ)*(pY-I*pX),
          ( pX + I*pY), 0 } );
      if( polState == -1 ) return N * Tensor({ pZ + pY*invPT2*(pP-pZ)*(pY+I*pX), 
          I*pZ - pX*invPT2*(pP-pZ)*(pY+I*pX),
          ( -pX + I*pY), 0 } );
      if( polState == 0 ) ERROR("Photon should does not have a rest frame, cannot be in m=0 state");
    }
  }

  if( m_props->twoSpin() == 1 && basis == "Weyl" ){
    Expression n       = fcn::sqrt( 2 * pP*(pP+pZ) );
    Expression fa      = fcn::sqrt( (pE + m)/(2*m) );
    Expression fb      = fcn::sqrt( (pE - m)/(2*m) );
    Expression aligned = make_cse( Abs(pP + pZ) < 10e-6 ) ;
    
    Expression xi10    = make_cse(Ternary( aligned,  1, (pP+pZ)/n ));
    Expression xi11    = make_cse(Ternary( aligned,  0,  z/n ));
    Expression xi00    = make_cse(Ternary( aligned,  0, -zb/n ));
    Expression xi01    = make_cse(Ternary( aligned,  1, (pP+pZ)/n ));


    if(id > 0 && polState ==  1 ) return Tensor({ fa*xi10,  fa*xi11,  fb*xi10,  fb*xi11 } );
    if(id > 0 && polState == -1 ) return Tensor({ fa*xi00,  fa*xi01, -fb*xi00, -fb*xi01 } );
    if(id < 0 && polState ==  1 ) return Tensor({ fb*xi00,  fb*xi01,  -fa*xi00,  -fa*xi01 } );
    if(id < 0 && polState == -1 ) return Tensor({ fb*xi10,  fb*xi11, -fa*xi01,  -fa*xi11 } );
  } 
  if ( m_props->twoSpin() == 1 && basis == "Dirac" )
  {
    Expression norm = make_cse( fcn::sqrt((pE+m)/(2*m) ));
    if(id > 0 && polState ==  1 ) return norm * Tensor({ 1        ,0          , pZ/(pE+m),  z/(pE+m)  });
    if(id > 0 && polState == -1 ) return norm * Tensor({ 0        ,1          , zb/(pE+m), -pZ/(pE+m) }); 
    if(id < 0 && polState == -1 ) return norm * Tensor({ pZ/(pE+m),  z/(pE+m) , 1        ,0           });
    if(id < 0 && polState ==  1 ) return norm * Tensor({ zb/(pE+m),-pZ/(pE+m) , 0        ,1           }); 
  }
  WARNING("Spin tensors not implemented for spin J = " << m_props->twoSpin() << " / 2 m = " << m_polState << " " << m_name ); 
  return Tensor( std::vector<double>( {1.} ), std::vector<size_t>( {0} ) );
}

bool Particle::checkExists() const
{
  bool success = true;
  if ( m_daughters.size() == 2 ) {
    success &= VertexFactory::isVertex( vertexName() );
    if ( !success ) {
      ERROR( uniqueString() );
      ERROR( "Spin configuration not found J = "
          << spin() << " L  =  " << m_orbital << " l' = " << m_spinConfigurationNumber
          << " s1 = " << m_daughters[0]->spin() << " s2 = " << m_daughters[1]->spin() );
    }
    success &= daughter( 0 )->checkExists() & daughter( 1 )->checkExists();
  }
  if ( success == false ) ERROR( uniqueString() << " is not described in IVertex" );
  return success;
}

std::pair<size_t, size_t> Particle::orbitalRange( const bool& conserveParity ) const
{
  if ( m_daughters.size() == 0 ) return {0, 0};
  if ( m_daughters.size() != 2 ) {
    ERROR( "L not well defined for D == " << m_daughters.size() );
    return {999, 998};
  }
  const int S  = m_props->twoSpin();
  const int s1 = daughter( 0 )->props()->twoSpin();
  const int s2 = daughter( 1 )->props()->twoSpin();

  int min                                  = std::abs( S - s1 - s2 );
  if ( std::abs( S + s1 - s2 ) < min ) min = std::abs( S + s1 - s2 );
  if ( std::abs( S - s1 + s2 ) < min ) min = std::abs( S - s1 + s2 );
  int max                                  = S + s1 + s2;

  min /= 2;
  max /= 2;
  DEBUG( "Range = " << min << " -> " << max << " conserving parity ? " << conserveParity << " J = " << S << " s1= " << s1 << " s2= " << s2 );
  if ( conserveParity == false ) return {min, max};
  
  int l = min;
  for ( ; l < max + 1; ++l ) if( conservesParity(l) ) break;
  if ( l == max + 1 ) return {999, 999};
  std::pair<size_t, size_t> lLimit = {l, l};
  l = max;
  for ( ; l != min - 1; --l )
    if ( conservesParity( l ) ) break;
  lLimit.second = l;
  return lLimit;
}

std::vector< std::pair<double,double> > Particle::spinOrbitCouplings( const bool& conserveParity ) const 
{   
  auto lRange = orbitalRange( conserveParity );
  std::vector< std::pair<double, double> > lsCouplings; 
  const double J  = m_props->twoSpin() / 2.0;
  const double j1 = daughter( 0 )->props()->twoSpin() /2.0;
  const double j2 = daughter( 1 )->props()->twoSpin() /2.0;
  for( double L = lRange.first; L <= lRange.second; ++L )
  {
    if( conserveParity && !conservesParity(L) ) continue; 
    for( double S = std::abs(j1-j2); S <= j1+j2 ; ++S ){
      double z2 = J*J - L*L - S*S;
      if( L != 0 && S != 0  && std::abs(z2/(2*L*S)) <= 1 ) lsCouplings.emplace_back(L,S);
      if( (L == 0 || S == 0) && z2 == 0 )                  lsCouplings.emplace_back(L,S);
    }
  }
  return lsCouplings; 
}

std::string Particle::texLabel( const bool& printHead, const bool& recurse ) const
{
  const std::string leftBrace     = "\\left[";
  const std::string rightBrace    = "\\right]";
  if( m_daughters.size() == 0 ) return m_props->label();
  
  std::string val                 = m_props->label();
  if( printHead && m_isHead ) val = m_props->label() + "\\rightarrow ";

  if( m_isHead || recurse ){
    if ( !m_isHead || m_orbital != m_minL ) val += leftBrace;
    for( auto& prod : m_daughters )         val += prod->texLabel(false,recurse) + " ";
    if ( !m_isHead || m_orbital != m_minL ) val += rightBrace;
    if ( m_orbital != m_minL )              val += "^{" + orbitalString() + "}";
  }
  return val;
}

void Particle::sortDaughters()
{
  std::stable_sort( m_daughters.begin(), m_daughters.end(), []( auto& A, auto& B ) { return *A < *B; } );
}

int Particle::conjugate( bool invertHead , bool reorder )
{
  int sgn = 1;
  DEBUG( "Conjugating flavour of " << m_name );
  if ( m_props == nullptr ) ERROR( "Particle properties not found for " << m_name );
  if ( ( !m_isHead || invertHead ) && m_props->hasDistinctAnti() ) {
    DEBUG( "Trying to conjugate: " << m_name );
    m_props = ParticlePropertiesList::get( -m_props->pdgID() );
    m_name  = m_props->name();
    m_parity = m_props->P();
  }
  sgn *= pow( -1, m_orbital );
  if ( m_parity == -1 
      && m_daughters.size() == 2 
      && std::abs(daughter(0)->props()->pdgID()) == std::abs(daughter(1)->props()->pdgID())
      && daughter(0)->props()->pdgID() != daughter(1)->props()->pdgID() ) {
    sgn *= -1;
  }
  for ( auto& d : m_daughters ) sgn *= d->conjugate( invertHead );
  if( reorder ) sortDaughters();
  std::string oldUniqueString = m_uniqueString;
  m_uniqueString              = makeUniqueString();
  DEBUG( "Conjugate : " << oldUniqueString << " to " << m_uniqueString << " sgn = " << sgn );
  return sgn;
}

QuarkContent Particle::quarks() const
{
  return m_name.find("NonRes") != std::string::npos ? daughterQuarks() : m_props->netQuarkContent();
}

QuarkContent Particle::daughterQuarks() const
{
  QuarkContent quarks;
  for ( auto& d : m_daughters ) quarks += d->quarks();
  return quarks;
}

Expression Particle::massSq() const { return isStable() ? mass() * mass() : make_cse( dot( P(), P() ) ); }
void Particle::addModifier( const std::string& mod )
{
  m_modifiers.push_back( mod );
  parseModifier( mod );
  m_uniqueString = makeUniqueString();
}
void Particle::setTop( bool state ) { m_isHead = state; }
void Particle::setIndex( const unsigned int& index, const bool& setOri )
{
  m_index                       = index;
  if ( setOri ) m_originalIndex = m_index;
}

void Particle::setOrbital( const unsigned int& orbital )
{
  m_orbital      = orbital;
  m_uniqueString = makeUniqueString();
}
void Particle::setLineshape( const std::string& lineshape )
{
  m_lineshape    = lineshape;
  m_uniqueString = makeUniqueString();
}
int Particle::finalStateParity() const
{
  int lpart = ( m_orbital % 2 == 0 ? 1 : -1 );
  for ( auto& d : m_daughters ) lpart *= d->parity();
  return lpart;
}

bool Particle::conservesParity( unsigned int L ) const
{
  return parity() == daughter( 0 )->parity() * daughter( 1 )->parity() * ( L % 2 == 0 ? 1 : -1 );
}

std::string Particle::topologicalString() const
{
  std::string topo = "";
  for ( auto& d : m_daughters ) topo += d->m_props->J();
  return topo;
}
const ParticleProperties* Particle::props() const { return m_props; }
bool Particle::isHead() const { return m_isHead; }
bool Particle::isWeakDecay() const { return quarks() == daughterQuarks(); }
bool Particle::isStateGood() const { return m_isStateGood; }
bool Particle::isStable() const { return m_daughters.size() == 0; }

bool Particle::isQuasiStable() const
{
  return props()->width() < ParticlePropertiesList::getMe()->quasiStableThreshold();
}
unsigned int Particle::orbital()  const { return m_orbital; }
int Particle::polState() const { return m_polState; }
std::string Particle::name() const { return m_name; }
std::string Particle::uniqueString() const { return m_uniqueString; }
std::string Particle::lineshape() const { return m_lineshape; }
unsigned int Particle::index() const { return m_index; }
unsigned int Particle::originalIndex() const { return m_originalIndex; }

int Particle::parity() const { return m_parity; }
double Particle::mass() const { return m_props->mass(); }
void Particle::setDaughter( const Particle& particle, const unsigned int& index )
{
  m_daughters[index] = std::make_shared<Particle>( particle );
  makeUniqueString();
}

std::vector<std::shared_ptr<Particle>> Particle::daughters() const { return m_daughters; }
bool Particle::operator<( const Particle& other )
{
  if ( spin() != other.spin() ) return spin() > other.spin();
  if ( !isStable() && other.isStable() ) return true;
  if ( isStable() && !other.isStable() ) return false;
  if ( mass() != other.mass() ) return mass() > other.mass();
  if ( std::abs( props()->pdgID() )  == std::abs( other.props()->pdgID() )
      && props()->pdgID() != other.props()->pdgID() ) return props()->pdgID() > other.props()->pdgID();
  if ( props()->charge() != other.props()->charge() ) return props()->charge() > other.props()->charge();
  if ( props()->pdgID() != other.props()->pdgID() ) return props()->pdgID() > other.props()->pdgID();
  return index() < other.index();
}

bool Particle::operator>(const Particle& other )
{
  return !( *this < other );
}

EventType Particle::eventType() const
{
  std::vector<std::string> names = {m_name};
  auto fs                        = getFinalStateParticles( false );
  for ( auto& p : fs ) names.push_back( p->name() );
  return EventType( names );
}

void Particle::setPolarisationState( const int& state )
{
  if( m_props->isFermion() && (state > m_props->twoSpin() || state == 0 ) ){
    WARNING("Invalid polarisation state for fermion ["<<name() << "] = " << state );
    m_polState = 999; 
  }
  else if ( m_props->isBoson() && state > m_props->twoSpin() / 2 ){
    WARNING("Invalid polarisation state for boson ["<< name() <<"] = " << state );
    m_polState = 999 ;
  }
  else {
    m_polState = state; 
  }
  m_uniqueString = makeUniqueString();
}

unsigned int Particle::matches( const Particle& other ) const 
{
  unsigned int rt=0;
  if( name() != other.name() ) return MatchState::None;
  if( (other.daughters().size() == 0 && daughters().size() != 0 ) || 
      (other.daughters().size() != 0 && daughters().size() == 0 ) ) return MatchState::PartialExpansion; 
  if( other.daughters().size() != daughters().size() ) return MatchState::None;
  for( unsigned int i = 0 ; i < daughters().size(); ++i )
  {
    auto daughter_match = daughters()[i]->matches( *other.daughters()[i] );
    if( daughter_match == 0 ) return MatchState::None;
    rt |= daughter_match;
  }
  if( m_orbital  != other.orbital() )  rt |= MatchState::DifferentOrbital;
  if( m_polState != other.polState() ) rt |= MatchState::DifferentPolarisation;
  if( rt == 0 ) rt = MatchState::Exact;
  if( rt & MatchState::Exact && rt != MatchState::Exact ) 
    rt &= ~MatchState::Exact;
  return rt;
}

std::string Particle::decayDescriptor() const { return m_uniqueString ; }

int Particle::quasiCP() const 
{
  int prod = m_props->C() == 0 ? 1 : m_props->C();
  prod *= ( m_orbital % 2 == 0 ? 1 : -1 );
  for( auto& d : m_daughters ) prod *= d->quasiCP() ;
  return prod; 
}

std::pair<bool,std::string> Particle::attribute(const std::string& key) const 
{
  std::map<std::string,std::string> m_rt;
  for( auto& modifier : m_modifiers ){
    auto tokens = split( modifier, '=');
    if( tokens.size() == 2 && tokens[0] == key ) return std::make_pair(true, tokens[1] );
  }
  return std::make_pair(false, "");
}
