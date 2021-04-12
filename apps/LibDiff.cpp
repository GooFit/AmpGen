#include "AmpGen/PhaseSpace.h"
#include "AmpGen/DynamicFCN.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

std::vector<std::string> cmd( const std::string& command ){
  std::vector<std::string> output;
  FILE* proc = popen(command.c_str(), "r");
  char buf[4096];
  while (!feof(proc) && fgets(buf, sizeof(buf), proc))  output.push_back(buf);
  pclose(proc);
  return output;
}

struct Counter {
  int pass {0};
  int total{0};
  Counter() = default;    
  Counter( const bool& pass ) : pass(pass), total(1) {}
  Counter( const double& a, const double& b, const std::string& name, const double& tolerance );
  Counter( const complex_t& a, const complex_t& b, const std::string& name, const double& tolerance );
  Counter(const std::vector<complex_t>& a, const std::vector<complex_t>& b, const std::string& name, const double& tolerance);
  void operator+=( const Counter& other )
  {
    this->pass  += other.pass;
    this->total += other.total;
  }

};

Counter::Counter(const double& a, const double& b, const std::string& name, const double& tolerance)
{
  double diff = std::abs(a-b);
  total = 1;
  if( diff > std::abs(tolerance) ){
    ERROR( name << " (" << a << " " << b << ") = " << diff << " > " << tolerance); 
    pass = 0;
  }
  else pass =1;
  DEBUG( name << " (" << a << " " << b << ") = " << diff  << " > " << tolerance); 
}

Counter::Counter(const complex_t& a, const complex_t& b, const std::string& name, const double& tolerance)
{
  double diff_re = std::abs(std::real(a)-std::real(b));
  double diff_im = std::abs(std::imag(a)-std::imag(b));
  total =1;
  if( diff_re > std::abs(tolerance)  || diff_im > std::abs(tolerance) ){
    ERROR( name << " (" << a << " " << b << ") = " << diff_re << ", " << diff_im << " > " << tolerance); 
    pass =0;
  }
  else pass = 1;
  DEBUG( name << " (" << a << " " << b << ") = " << diff_re << ", " << diff_im << " < " << tolerance); 
}

Counter::Counter(const std::vector<complex_t>& a, const std::vector<complex_t>& b, const std::string& name, const double& tolerance)
{
  total = a.size();
  for( size_t i = 0 ; i < a.size(); ++i ){
    double diff_re = std::abs(std::real(a[i])-std::real(b[i]));
    double diff_im = std::abs(std::imag(a[i])-std::imag(b[i]));
    if( diff_re > std::abs(tolerance)  || diff_im > std::abs(tolerance) ){
      ERROR( name << " (" << a[i] << " " << b[i] << ") = " << diff_re << ", " << diff_im << " > " << tolerance); 
    }
    else pass++;
    DEBUG( name << " (" << a[i] << " " << b[i] << ") = " << diff_re << ", " << diff_im << " > " << tolerance); 
  }
}

int main( int argc, char** argv )
{
  std::string modelName = argv[1];
  OptionsParser::getMe()->setQuiet();
  OptionsParser::setArgs( argc, argv );
  auto t = EventType(NamedParameter<std::string>("EventType","").getVector()); 
  PhaseSpace phsp(t);
  auto event = phsp.makeEvent();
  auto event2 = phsp.makeEvent();
  std::string lib    = NamedParameter<std::string>("Lib","");
  std::string refLib = NamedParameter<std::string>("RefLib","");
  std::string type   = NamedParameter<std::string>("Type","CoherentSum");
  Counter total;

  auto ftable = cmd("nm " + refLib + "| grep __wParams");
  if( type == "CoherentSum"){
    auto f1  = DynamicFCN<complex_t(const double*, const int&)>(lib, "AMP");
    auto f2  = DynamicFCN<complex_t(const double*, const int&)>(refLib, "AMP");
    total += Counter(f1(event,+1), f2(event,+1), modelName + " fcn(x)", 10e-6 );
    total += Counter(f1(event,-1), f2(event,-1), modelName + " fcn(Px)", 10e-6);
    for( auto& line : ftable )
    {
      auto i  = trim( split(line,' ')[2] );
      auto f  = DynamicFCN<complex_t(const double*)>(lib, i );
      auto g  = DynamicFCN<complex_t(const double*)>(refLib, i );
      if( !f.isLinked() || !g.isLinked() ) total += false ;
      else total += Counter( f(event), g(event), i, 10e-6);
    }
  }
  if( type == "PolarisedSum"){
    if( std::get<0>( t.dim() ) > 1 )
    {
      auto f1  = DynamicFCN<double(const double*, const int&, const double&, const double&, const double&)>(lib   , "FCN_extPol");
      auto f2  = DynamicFCN<double(const double*, const int&, const double&, const double&, const double&)>(refLib, "FCN_extPol");
      total += Counter(f1(event,+1,0,0,0), f2(event,+1,0,0,0), modelName + " fcn(x)" , 1e-6);
      total += Counter(f1(event,-1,0,0,0), f2(event,-1,0,0,0), modelName + " fcn(Px)", 1e-6);
    }
    for( auto& line : ftable ){
      auto i       = trim( split(line,' ')[2] );
      auto a1  = DynamicFCN<std::vector<complex_t>(const double*)>(lib   , i );
      auto a2  = DynamicFCN<std::vector<complex_t>(const double*)>(refLib, i );
      if( ! a1.isLinked() || ! a2.isLinked() ) total += false; 
      else {
        auto pass = Counter(a1(event), a2(event), i, 10e-6);
        total += pass;
      }
    }
  }
  if(total.pass == total.total) 
    INFO("Library: " << modelName << " matches [passed = " << total.pass <<  " / " << total.total << "]" );
  else 
    ERROR("Library: " << modelName << " does not match [passed = " << total.pass <<  " / " << total.total << "]" );
}
