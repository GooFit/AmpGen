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
using namespace std::complex_literals; 

Particle::Particle() = default; 

Particle::Particle( const std::string& name, const Particle& p1, const Particle& p2 ) :
   m_name{name}
 , m_props{ParticlePropertiesList::get(name)}
 , m_daughters{ std::make_shared<Particle>(p1), std::make_shared<Particle>(p2) }
{
  pdgLookup();
  sortDaughters();
  for ( auto& d : m_daughters ) d->setParent( this );
  m_uniqueString = makeUniqueString();
}

Particle::Particle( const int& pdgID, const Particle& p1, const Particle& p2 ) :
    m_props{ParticlePropertiesList::get(pdgID)}
  , m_daughters{std::make_shared<Particle>(p1) , std::make_shared<Particle>(p2) }
{ 
  if ( m_props != nullptr ) m_name = m_props->name();
  pdgLookup();
  sortDaughters();
  for ( auto& d : m_daughters ) d->setParent( this );
  m_uniqueString = makeUniqueString();
}

Particle::Particle( const std::string& name, const unsigned int& index ) :
   m_props{ParticlePropertiesList::get( name )}
 , m_index{index}
{
  pdgLookup();
  m_uniqueString = makeUniqueString();
}

Particle::Particle( const std::string& decayString, const std::vector<std::string>& finalStates,
    const bool& orderDaughters )
{
  auto items = getItems( decayString );
  if ( items.size() == 0 ) {
    ERROR( "Failed to parse string" << decayString );
    m_isStateGood = false;
    return;
  }
  m_modifiers = getItems( items[0], {"[", "]"}, ";" );
  m_name      = m_modifiers[0];
  m_modifiers.erase( m_modifiers.begin() );
  for( auto& modifier : m_modifiers ) parseModifier( modifier );
  for( size_t it = 1; it < items.size(); ++it ){
    m_daughters.push_back( std::make_shared<Particle>( items[it] ) );
    (*m_daughters.rbegin())->setParent(this);
  }
  m_props = ParticlePropertiesList::get( m_name );
  pdgLookup();
  auto fs = getFinalStateParticles( false );
  DEBUG( "Building particle from decay string = " << decayString << " name = " << m_name << " "<< m_props->name() );
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
    m_isStateGood = std::all_of( hasUsedFinalState.begin(), hasUsedFinalState.end(),[](const auto& b){return b;} );
  }
  m_uniqueString = makeUniqueString();
}


bool Particle::isValidDecayDescriptor( const std::string& decayDescriptor )
{
  size_t firstOpen = decayDescriptor.find("{");
  size_t open  = std::count(decayDescriptor.begin(), decayDescriptor.end(), '{'); 
  size_t close = std::count(decayDescriptor.begin(), decayDescriptor.end(), '}'); 
  if( open == 0 || open != close || firstOpen == std::string::npos) return false; 
  std::string firstState = decayDescriptor.substr(0, firstOpen);
  auto firstSquare = firstState.find("[");
  if( firstSquare == std::string::npos ) return ParticleProperties::get( firstState, true ) != nullptr;
  return ParticleProperties::get( firstState.substr(0, firstSquare), true ) != nullptr; 
}

void Particle::parseModifier( const std::string& mod )
{
  if ( Lineshape::Factory::isLineshape(mod) ) m_lineshape = mod;
  else if( mod.size() == 1 )
  {
    DEBUG( "Modifier = " << mod );
    if ( mod == "S" ) m_orbital = 0;
    else if ( mod == "P" ) m_orbital = 1;
    else if ( mod == "D" ) m_orbital = 2;
    else if ( mod == "F" ) m_orbital = 3;
    else if ( mod == "G" ) m_orbital = 4;
    else if ( mod == "H" ) m_orbital = 5;
    else if ( mod == "I" ) m_orbital = 6;
    else { 
      bool status=true;
      auto tmp = lexical_cast<unsigned int>(mod,status); // this is if it is a number /// 
      if( status ) m_spinConfigurationNumber = tmp;
      else {
        m_spinConfigurationNumber = 10 + int(mod[0]);
      }
    }
  }
  else if ( mod.size() == 2 )
  {
    parseModifier( mod.substr(0,1) );
    parseModifier( mod.substr(1,1) );
  }
}

double Particle::spin() const { return double( m_props->twoSpin() / 2. ) ; }
double Particle::S() const { return m_spinConfigurationNumber ; }

void Particle::pdgLookup()
{
  if ( m_props == nullptr ) {
    m_isStateGood = false;
    return;
  }
  if ( m_lineshape == "BW" || m_usesDefaultLineshape ) {
    if ( m_name.find("NonRes") != std::string::npos || m_props->width() < ParticlePropertiesList::getMe()->quasiStableThreshold() ) m_lineshape = "FormFactor";
    if ( m_props->isPhoton() ) m_lineshape = "Photon";
    m_usesDefaultLineshape = true;
  } 
  else m_usesDefaultLineshape = false;
  m_parity                 = m_props->P();
  // if ( !isdigit( m_props->J()[0] ) ) ERROR( "Spin not recognised! : " << m_name << " J = " << m_props->J() );
  if( m_defaultModifier != "" && m_lineshape.find(".") == std::string::npos ){
    m_lineshape = m_lineshape + "." + m_defaultModifier.getVal();
  }
  bool isStrong = quarks() == daughterQuarks();
  if( abs(m_props->pdgID()) == 24 || abs(m_props->pdgID()) == 23 ) isStrong = false; 
  if ( m_name.find( "NonRes" ) != std::string::npos ) isStrong = true;
  m_minL = m_daughters.size() == 2 ? orbitalRange( isStrong ).first : 0;
  if ( m_daughters.size() == 2 ) {
    DEBUG( "IsStrong ? " << isStrong << " orbital =  " << m_orbital << " minL   =  " << m_minL
        << " name   =  " << m_name << " d0     =  " << m_daughters[0]->name()
        << " d1     =  " << m_daughters[1]->name() );
  }
  if ( m_orbital == 0 ) m_orbital = m_minL; // define in ground state unless specified
  if( m_daughters.size() != 0 ){
    DEBUG( m_name << " is decaying via " << ( isStrong ? "strong" : "electroweak" ) << " interactions; P = " << props()->P() << "l = " << m_orbital );
  }
  int charge = 0; 
  for ( auto& d : m_daughters ){
    d->setParent(this);
    charge += d->props()->charge();
  }
  if( m_minL == 999 ) ERROR("Decay: " << m_name << " does not appear to have an allowed spin-orbit configuration");
  if( m_daughters.size() != 0 && m_props->charge() != charge ) ERROR("Decay: " << m_name << " does not conserve (electric) charge");
}

Tensor Particle::P() const
{
  Tensor momentum( std::vector<double>( {0., 0., 0., 0.} ), Tensor::dim(4) );
  if ( isStable() ) {
    if ( m_index != 999 ) {
      const std::string index = std::to_string( m_index );
      Tensor rt( std::vector<Expression>( {
            Parameter(index + "_Px"), 
            Parameter(index + "_Py"),
            Parameter(index + "_Pz"), 
            Parameter(index + "_E") }) , Tensor::dim(4) );
//      rt[3] = fcn::sqrt( mass()*mass() + rt[0]*rt[0] + rt[1]*rt[1] + rt[2]*rt[2] ) ;
      return rt;
    } else ERROR( "Stable particle " << m_index << "is unindexed!" );
  } 
  else {
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

std::shared_ptr<Particle> Particle::daughter( const std::string& name, const int& maxDepth ) const 
{
  if( maxDepth == -1 )
  {
    auto fs = getFinalStateParticles();
    for( auto& f : fs ) if( f->name() == name ) return f;
    ERROR("Particle: " << name << " not found amongst decay products!");
  }
  else {
    for( auto& d : m_daughters ) if( d->name() == name ) return d;
    ERROR("Particle: " << name << " not found amongst decay products!");
  }
  return std::shared_ptr<Particle>();
}


std::string Particle::orbitalString() const
{
  constexpr std::array<char, 7> orbitals = {'S','P','D','F','G','H','I'};
  std::string rt = std::string(1, orbitals[m_orbital] ); 
  if( m_spinConfigurationNumber != 0 ){ 
    if( m_spinConfigurationNumber < 10 ) rt += std::to_string( m_spinConfigurationNumber  );
    else rt += char(m_spinConfigurationNumber-10);
  }
  return rt; 
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
    m_uniqueString = m_props->name() + modifier + "{" + vectorToString( m_daughters, ",", [](auto& d){ return d->makeUniqueString() ;} ) +"}";
  } else {
    m_uniqueString = m_props->name() + modifier ;
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
    if( qTree.daughters().size() == 0 || (qTree.isQuasiStable() && m_daughters.size() != 1 ) )      
      returnVal.addDaughter( std::make_shared<Particle>(qTree) );
    else 
      for( auto& it : qTree.daughters() ) returnVal.addDaughter(it); 
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
  Expression s = massSq();
  Expression prop = 1; 
  if ( m_daughters.size() == 2 ) 
    prop = Lineshape::Factory::get(m_lineshape, s, daughter(0)->massSq(), daughter(1)->massSq(), m_name, m_orbital, db);
  else if ( m_daughters.size() >= 3 )
    prop = Lineshape::Factory::get(m_lineshape == "BW" ? "SBW" : m_lineshape, *this, db );
  else if ( m_daughters.size() == 1 && m_lineshape != "BW" && m_lineshape != "FormFactor" )
  { 
    prop = Lineshape::Factory::get(m_lineshape, *this, db );
  } 
  total = total * make_cse(prop);
  for(auto& d : m_daughters) total = total*make_cse(d->propagator(db));
  if(db != nullptr) db->emplace_back("A("+uniqueString()+")", total);
  return total;
}

std::vector<std::vector<size_t>> Particle::identicalDaughterOrderings() const
{
  std::vector<std::vector<size_t>> orderings;
  auto finalStateParticles = getFinalStateParticles();
  std::vector<std::string> zero_perm;
  std::transform( finalStateParticles.begin(), finalStateParticles.end(), std::back_inserter(zero_perm), 
      [](auto& p ){ return p->name() ; } );
  std::vector<size_t> indices( finalStateParticles.size() );
  std::iota(indices.begin(), indices.end(), 0);
  do {
    bool isSameAsP0 = true;
    for ( size_t i = 0; i < indices.size(); ++i )
      if ( zero_perm[i] != finalStateParticles[indices[i]]->name() ) isSameAsP0 = false;
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
void Particle::addDaughter( const std::shared_ptr<Particle>& particle ) 
{ 
  m_daughters.push_back( particle ); 
  m_uniqueString = makeUniqueString();
}

Tensor Particle::transitionMatrix( DebugSymbols* db  )
{
  auto particles = getFinalStateParticles();
  return spinTensor();
}

Expression Particle::getExpression( DebugSymbols* db, const std::vector<int>& state)
{
  if( state.size() !=0 ) setPolarisationState( state );
  if( db != nullptr && !isStable() ) 
    db->emplace_back( uniqueString() , Parameter( "NULL", 0, true ) );
  Expression total = 0;
  Tensor::Index a;
  auto finalStateParticles  = getFinalStateParticles();
  auto orderings            = identicalDaughterOrderings();
  bool doSymmetrisation = !hasModifier("NoSym");
  bool sumAmplitudes    = !hasModifier("Inco");
  bool includeSpin      = !hasModifier("BgSpin0");
  std::vector<int> exchangeParities;
  std::transform( finalStateParticles.begin(), finalStateParticles.end(), std::back_inserter(exchangeParities), 
      [](const auto& p){ return p->props()->isFermion() ? -1 : 1; } );
  for(const auto& ordering : orderings){
    auto exchangeParity = minSwaps( ordering, exchangeParities );   
    setOrdering( ordering );
    Expression spinFactor = 1; 
    if( includeSpin && m_spinFormalism == spinFormalism::Covariant ){
      Tensor st = spinTensor(db);
      if( m_props->twoSpin() == 0 ) spinFactor = st[0];
      if( m_props->twoSpin() == 1 ){
        Tensor is = Bar( externalSpinTensor(m_polState) );
        ADD_DEBUG_TENSOR( externalSpinTensor(m_polState), db );
        ADD_DEBUG_TENSOR( st, db );
        if( st.size() != 4 ){ 
          ERROR("Spin tensor is the wrong rank = " << st.dimString() ); 
          spinFactor = 1; 
        }
        else spinFactor = is(a) * st(a) ;
      }
      if( m_props->twoSpin() == 2 ){
        Tensor is  = externalSpinTensor(m_polState).conjugate();
        ADD_DEBUG_TENSOR(is, db );
        spinFactor = dot(is,st);
      }
    }
    if ( includeSpin && m_spinFormalism == spinFormalism::Canonical ){
      spinFactor = helicityAmplitude(*this, TransformSequence(), m_props->isBoson() ? polState() : double(polState())/2.0, db);
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
    Expression ls = make_cse(propagator(db)) * spinFactor; 
    if(sumAmplitudes) total = total + exchangeParity.second * ls;
    else              total = total + fcn::conj(ls) * ls;
    if(!doSymmetrisation) break;
  }
  ADD_DEBUG( total, db );
  double nPermutations = doSymmetrisation ? orderings.size() : 1;
  if ( sumAmplitudes )
  {
    if ( is<Constant>(total) ){
      WARNING("Amplitude is just a constant: " << total << " may cause problems for compiler, making a little bit complex" );
      total += 1i * 0.00001;
    }
    return total / fcn::sqrt( nPermutations );
  }
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
    auto vname = m_props->spinName() + "_" + m_daughters[0]->m_props->spinName() + m_daughters[1]->m_props->spinName() + "_" + orbitalString();
    Tensor value = Vertex::Factory::getSpinFactor( P(), Q(), 
        daughter(0)->spinTensor(db),
        daughter(1)->spinTensor(db), vname, db );
    DEBUG( "Returning spin tensor" );
    return value;
  } else if ( m_daughters.size() == 3 ) {
    return Vertex::Factory::getSpinFactorNBody( {
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
    return Tensor( std::vector<double>( {1.} ), std::vector<unsigned>( {1} ) );
  if( m_spinBasis == spinBasis::Dirac && ( m_props->isPhoton() || m_props->isNeutrino() ) )
  {
    WARNING("Use the Weyl (helicity) basis for calculations involving photons / neutrinos, as they don't have a rest frame. This will result in ill-defined amplitudes.");
  }
  
  Tensor p        = P();
  Expression pX   = p.get(0);
  Expression pY   = p.get(1);
  Expression pZ   = p.get(2);
  Expression pE   = p.get(3);
  Expression pP   = fcn::sqrt( pX*pX + pY*pY + pZ*pZ );
  Expression m    = mass();
  Expression z    = pX + 1i*pY; 
  Expression zb   = pX - 1i*pY;
  auto id = props()->pdgID(); 
  if ( m_props->twoSpin() == 2 )
  {
    if( m_spinBasis == spinBasis::Weyl ) 
    {
      Expression N = 1./(sqrt(2)); 
      Expression invPT2 = make_cse(Ternary( pX*pX + pY*pY > 1e-6, 1./(pX*pX+pY*pY), 0 ));
      Expression invP   = make_cse(Ternary( pP            > 1e-6, 1./pP,0));
      Expression pZoverP = make_cse(Ternary( pP > 1e-6, pZ/pP, 1) );
      Expression f = (pP-pZ)*invPT2;
      if(polState ==  1) return N * Tensor({-pZoverP + 1i*z *pY*f*invP, -1i*pZoverP - 1i*z *pX*f*invP, z*invP           , 0.  }); 
      if(polState == -1) return N * Tensor({ pZoverP + 1i*zb*pY*f*invP, -1i*pZoverP - 1i*zb*pX*f*invP,-zb*invP          , 0.  });
      if(polState ==  0) return     Tensor({pX*pE*invP/m    , pY*pE*invP/m       , pZoverP*pE/m, pP/m}); 
    }
    if( m_spinBasis == spinBasis::Dirac ){
      Expression N = make_cse(1./(m*(pE + m)));
      if( polState ==  1 ) return -Tensor({1.+ z *pX*N,  1i +  z*pY*N,  z*pZ*N    ,  z*m })/sqrt(2);
      if( polState == -1 ) return  Tensor({1.+ zb*pX*N, -1i + zb*pY*N, zb*pZ*N    , zb*m })/sqrt(2);
      if( polState ==  0 ) return  Tensor({pX*pZ*N    ,       pY*pZ*N, 1 + pZ*pZ*N, pZ/m });
    } 
  }
  if( m_props->twoSpin() == 1 )
  {
    if( m_spinBasis == spinBasis::Weyl )
    {
      std::array<Tensor,2> xi;
      Expression n       = fcn::sqrt( 2 * pP*(pP+pZ) );

      xi[0] = Tensor( {-zb/n , (pP+pZ)/n});
      xi[1] = Tensor( {(pP+pZ)/n, z/n });
      
      Expression fa = m_props->isNeutrino() ? polState * fcn::sqrt(pE) : polState * fcn::sqrt( pE/m-  1 );
      Expression fb = m_props->isNeutrino() ? fcn::sqrt(pE)            : fcn::sqrt( pE/m + 1 );
      int ind = (id>0 ? polState : -polState) == -1 ? 0 : 1;
     
      return id > 0 ?  Tensor({ - fa * xi[ind][0], - fa * xi[ind][1], fb * xi[ind][0], fb * xi[ind][1] })
        :  -polState * Tensor({   fb * xi[ind][0],   fb * xi[ind][1], fa * xi[ind][0], fa * xi[ind][1] });
    } 
    if ( m_spinBasis == spinBasis::Dirac )
    {
      Expression N = make_cse( fcn::sqrt((pE+m)/(2*m) ));
      if(id > 0 && polState ==  1 ) return N * Tensor({ 1        ,0          , pZ/(pE+m),  z/(pE+m)  });
      if(id > 0 && polState == -1 ) return N * Tensor({ 0        ,1          , zb/(pE+m), -pZ/(pE+m) }); 
      if(id < 0 && polState == -1 ) return N * Tensor({ pZ/(pE+m),  z/(pE+m) , 1        ,0           });
      if(id < 0 && polState ==  1 ) return N * Tensor({ zb/(pE+m),-pZ/(pE+m) , 0        ,1           }); 
    }
  }
  std::string js = m_props->isBoson() ? std::to_string(m_props->twoSpin()/2) : std::to_string(m_props->twoSpin()) +"/2";
  WARNING("Spin tensors not implemented for spin J = " << js << "; m = " << m_polState << " " << m_name ); 
  return Tensor( std::vector<double>( {1.} ), Tensor::dim(0) );
}

std::pair<size_t, size_t> Particle::orbitalRange( const bool& conserveParity ) const
{
  if( m_daughters.size() == 0 ) return {0, 0};
  if( m_daughters.size() == 1 ) return {0, 0};
  if( m_daughters.size() != 2 ) {
    ERROR( "L not well defined for nDaughters == " << m_daughters.size() );
    return {999, 998};
  }
  const int S  = m_props->twoSpin();
  const int s1 = daughter(0)->props()->twoSpin();
  const int s2 = daughter(1)->props()->twoSpin();
  int min = std::abs( S - s1 - s2 );
  min     = std::min(min, std::abs( S + s1 - s2 ));
  min     = std::min(min, std::abs( S - s1 + s2 ));
  int max = S + s1 + s2;
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

std::vector<std::pair<double,double>> Particle::spinOrbitCouplings( const bool& conserveParity ) const 
{   
  auto lRange = orbitalRange( conserveParity );
  std::vector<std::pair<double, double>> lsCouplings; 
  const double J  = m_props->twoSpin() / 2.0;
  const double j1 = daughter(0)->props()->twoSpin() /2.0;
  const double j2 = daughter(1)->props()->twoSpin() /2.0;
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
  std::string val                 = ""; // m_props->label();
  if( !isHead() )             val  = m_props->label();
  if( printHead && isHead() ) val = m_props->label() + "\\rightarrow ";
  
  if( isHead() || recurse ){
    if ( !isHead() || m_orbital != m_minL ) val += leftBrace;
    for( auto& prod : m_daughters )         val += prod->texLabel(false,recurse) + " ";
    if ( !isHead() || m_orbital != m_minL ) val += rightBrace;
    if ( m_orbital != m_minL )              val += "^{" + orbitalString() + "}";
  }
  return val;
}

void Particle::sortDaughters()
{
  bool isAntiParticle = m_props->pdgID() < 0;
  if( isAntiParticle ) conjThis();
  std::stable_sort( m_daughters.begin(), m_daughters.end(), [](const auto& A, const auto& B) { return *A < *B; } );
  if( isAntiParticle ) conjThis(); 
  m_uniqueString = makeUniqueString();
}

void Particle::conjThis()
{
  if( m_props->hasDistinctAnti() ) 
    m_props = ParticlePropertiesList::get( -m_props->pdgID(), false);
  for(unsigned i = 0 ; i != m_daughters.size(); ++i) m_daughters[i]->conjThis(); 
}

Particle Particle::conj( bool invertHead , bool reorder )
{
  Particle cp(*this);
  cp.clearDecayProducts();
  const ParticleProperties* p = ( !isHead() || invertHead ) && m_props->hasDistinctAnti() ? ParticlePropertiesList::get( -m_props->pdgID() ) : m_props;
  cp.setName( p->name() );
  for( auto& d : m_daughters ) cp.addDaughter( std::make_shared<Particle>(d->conj(invertHead,reorder)) );
  if( reorder ) cp.sortDaughters();
  cp.pdgLookup();
  return cp;
}

QuarkContent Particle::quarks() const
{
  return m_name.find("NonRes") != std::string::npos ? daughterQuarks() : m_props->quarkContent();
}

QuarkContent Particle::daughterQuarks() const
{
  QuarkContent quarks;
  for ( auto& d : m_daughters ) quarks += d->quarks();
  return quarks;
}

Expression Particle::massSq() const { return isStable() ? mass() * mass() : make_cse( dot( P(), P() ) ); }

void Particle::clearDecayProducts()
{
  m_daughters.clear();
}

void Particle::addModifier( const std::string& mod )
{
  m_modifiers.push_back( mod );
  parseModifier( mod );
  m_uniqueString = makeUniqueString();
}
void Particle::setName(const std::string& name)
{ 
  m_name   = name;
  m_props  = ParticlePropertiesList::get(name);
  m_parity = m_props->P();
}
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
  if( m_daughters.size() == 0) return 1; 
  if( m_daughters.size() == 1) return m_daughters[0]->finalStateParity(); /// forward through quasiparticles 
  if( m_daughters.size() == 2){
    int f = spin() - ( m_orbital + daughter(0)->spin() + daughter(1)->spin() );
    int lpart = (f % 2 == 0 ? 1 : -1);
    for( auto& d : m_daughters ) lpart *= d->finalStateParity();
    return lpart;
  }
  WARNING("> 2 body vertices may require special considerations when conjugating, returning (-1)^J");
  return pow(-1, spin() ); 
}

bool Particle::conservesParity( unsigned int L ) const
{
  return parity() == daughter(0)->parity() * daughter(1)->parity() * ( L % 2 == 0 ? 1 : -1 );
}

std::string Particle::topologicalString() const
{
  return std::accumulate( m_daughters.begin(), m_daughters.end(), std::string(""), [](auto& s, auto& d){ return s+d->props()->J() ; } );
}

const ParticleProperties* Particle::props() const { return m_props; }
bool Particle::isHead() const { return m_parent == nullptr; }
bool Particle::isWeakDecay() const { return quarks() != daughterQuarks(); }
bool Particle::isStateGood() const { return m_isStateGood; }
bool Particle::isStable() const { return m_daughters.size() == 0; }

bool Particle::isQuasiStable() const
{
  return props()->width() < ParticlePropertiesList::getMe()->quasiStableThreshold() && name() != "gamma0";
}
unsigned Particle::L()  const { return m_orbital; }
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

void Particle::setParent(const Particle* ptr )
{
  m_parent = ptr;
}

std::vector<std::shared_ptr<Particle>> Particle::daughters() const { return m_daughters; }

bool Particle::operator<( const Particle& other )
{
  if ( spin()      !=  other.spin()      ) return spin() > other.spin();
  if ( !isStable() &&  other.isStable()  ) return true;
  if ( isStable()  && !other.isStable()  ) return false;
  if ( mass()      !=  other.mass()      ) return mass() > other.mass();
  if ( std::abs(props()->pdgID())  == std::abs(other.props()->pdgID() )
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
  std::transform(fs.begin(), fs.end(), std::back_inserter(names), [](const auto& p){ return p->name(); } );
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

void Particle::setPolarisationState( const std::vector<int>& state )
{
  setPolarisationState( state[0] );
  auto fs = getFinalStateParticles();
  for(size_t i = 0 ; i < fs.size(); ++i ) fs[i]->setPolarisationState( state[i+1] );
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
  if( m_orbital  != other.L() )        rt |= MatchState::DifferentOrbital;
  if( m_polState != other.polState() ) rt |= MatchState::DifferentPolarisation;
  if( rt == 0 ) rt = MatchState::Exact;
  if( rt & MatchState::Exact && rt != MatchState::Exact ) 
    rt &= ~MatchState::Exact;
  return rt;
}

std::string Particle::decayDescriptor() const { return m_uniqueString ; }

int Particle::CP() const 
{
  return C() * finalStateParity();
}

int Particle::C() const 
{
  if( m_daughters.size() == 1 ) return m_daughters[0]->C();
  int prod = 1;
  for( auto& d : m_daughters ) prod *= ( d->props()->C() == 0 ? 1 : d->props()->C() ) * d->C();
  return prod; 
}

stdx::optional<std::string> Particle::attribute(const std::string& key) const 
{
  for( auto& modifier : m_modifiers ){
    auto tokens = split( modifier, '=');
    if( tokens.size() == 2 && tokens[0] == key ) return stdx::optional<std::string>(tokens[1]);
  }
  return stdx::nullopt;
}

std::ostream& AmpGen::operator<<( std::ostream& os, const Particle& particle ){
  return os << particle.decayDescriptor();
}

namespace AmpGen { 
  complete_enum( spinFormalism, Covariant, Canonical )
  complete_enum( spinBasis    , Dirac    , Weyl ) 
}

