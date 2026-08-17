#ifndef AMPGEN_OPTIONSPARSER_H
#define AMPGEN_OPTIONSPARSER_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <functional>
#include "AmpGen/MsgService.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen {

  class ConfigurableBase;

  class OptionsParser {
  public:
    typedef std::map<std::string, std::vector<std::string>>::const_iterator const_iterator;
    typedef std::map<std::string, std::vector<std::string>>::iterator iterator;

    static OptionsParser *getMe();
    static bool printHelp();
    static void setArgs(int argc, char **argv, const std::string &description = "");
    static void setArg(const std::string &arg);
    void setQuiet();
    void addArg(const std::string &arg);
    void setCommandLineArgs(int argc, char **argv, const std::string &description = "");
    void import(const std::string &fName);
    iterator find(const std::string &name);
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    std::vector<std::vector<std::string>> getInputOrdered() const;

    const std::vector<const ConfigurableBase *> &configurables() const { return m_configurables; }
    template <typename T> T *registerClass(T *object) {
      DEBUG("Registering -> " << type_string<T>());
      addToConfigurables(object);
      return object;
    }
    void print() const;
    void addToConfigurables(const ConfigurableBase *);

  private:
    std::vector<std::string> m_orderedKeys;
    std::map<std::string, std::vector<std::string>> m_parsedLines;
    std::map<std::string, std::function<void(std::vector<std::string>)>> m_keywords;
    std::vector<const ConfigurableBase *> m_configurables;

    bool m_printHelp = {false};
    bool m_quiet = {false};
    static OptionsParser *gOptionsParser;

    OptionsParser();
    bool ignoreThisLine(const std::string &line);
    void readStream(std::istream &is);
    std::vector<std::string> makeParsedStrings(const std::string &line, int &braceDepth) const;
    void addArg(const std::vector<std::string> &tokens);
  };
} // namespace AmpGen

#define REGISTER_CONFIGURABLE(CLASS_NAME) template <> CLASS_NAME *Configurable<CLASS_NAME>::gImpl = OptionsParser::getMe()->registerClass(new CLASS_NAME())

#endif
