#ifndef AMPGEN_PROGRESSBAR_H
#define AMPGEN_PROGRESSBAR_H

#include <string>
namespace AmpGen {
  class ProgressBar {
    public:
      ProgressBar(const size_t& width, const std::string& context);
      ~ProgressBar();
      void print(const double& percentage, const std::string& message="");
      void finish();
    private:
      size_t      m_width;
      int         m_lastPercent; 
      std::string m_context;
      std::string m_lastMessage = {""}; 
      bool        m_finished = {false};
  };
}

#endif
