#include "AmpGen/ProgressBar.h"
#include "AmpGen/MsgService.h"

#include <iostream>
#include <algorithm>
#include <iterator>
using namespace AmpGen;

ProgressBar::ProgressBar(const size_t& width, const std::string& context)
  : m_width(width),
  m_lastPercent(-1),
  m_context(context) {}

ProgressBar::~ProgressBar(){
  if( !m_finished ) finish();
}

void ProgressBar::print(const double& percentage, const std::string& message)
{
  int lpad = int(percentage * m_width);
  int val  = int(percentage * 100);
  if( val == m_lastPercent ) return; 
  m_lastPercent = val;
  std::cout << "\r\033[2;34m" << std::left << std::setw( FCNNAMELENGTH ) << m_context << "  INFO         "       << "\033[0m";
  std::cout << "Completed: "  << std::right << std::setw(3) << val << "% " << "[";
  std::fill_n(std::ostream_iterator<char>(std::cout), lpad, '|');
  std::fill_n(std::ostream_iterator<char>(std::cout), m_width-lpad, ' ');
  std::cout << "]";
  std::cout << message;
  m_lastMessage = message; 
  fflush (stdout);
}

void ProgressBar::finish(){
  print(1,m_lastMessage);
  std::cout << std::endl;
  m_finished = true; 
}
