#ifndef AMPGEN_PROPERTY_H
#define AMPGEN_PROPERTY_H
#include <stddef.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <map>
#include <cstring>
#include <string_view>

#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen
{
  /** Properties are configurable characteristics of objects that will at some point inherit from a base class
   *
   */
  template <typename value_t> class Property
  {
  protected:
    std::string m_name;
    std::string m_helpString;
    value_t m_value;

  public:
    // WARNING: NamedParameter style constructor will be depreciated at some point to allow for stricter typing
    Property(const std::string &name, const value_t &def = value_t(), const std::string_view &helpString = "")
        : m_name(name), m_helpString(helpString), m_value(def)
    {
      setFromOptionsParser();
      if(OptionsParser::printHelp())
        help(def);
      DEBUG(*this);
    }
    Property(void * /*parent*/, const std::string &name, const value_t &def = value_t(), const std::string_view &helpString = "")
        : m_name(name), m_helpString(helpString), m_value(def)
    {
      setFromOptionsParser();
      if(OptionsParser::printHelp())
        help(def);
      DEBUG(*this);
    }
    void set(const value_t &val){ m_value = val; } 
    template <typename T> bool operator==(const T &other) const { return m_value == other; }
    template <typename T> bool operator!=(const T &other) const { return m_value != other; }
    operator const value_t&() const { return m_value; }
    operator       value_t&()       { return m_value; }
    const value_t &value() const { return m_value; }
    const std::string &name() const { return m_name; }
    template <typename T> friend std::ostream &operator<<(std::ostream &os, const Property<T> &np);
    bool setFromStrings(const std::vector<std::string> &vsl)
    {
      bool status = true;
      if constexpr(isVector<value_t>::value)
        {
          m_value.resize(vsl.size());
          for(unsigned int i = 0; i < vsl.size(); i++)
            {
              m_value[i] = lexical_cast<typename value_t::value_type>(vsl[i], status);
              if(status == false)
                {
                  ERROR("Failed to parse token: " << vsl[i] << " for parameter: " << m_name);
                  return false;
                }
            }
        }
      else if constexpr(isTuple<value_t>::value)
        {
          for_each_with_counter(m_value, [this, vsl](auto &f, unsigned i) {
            bool status = true;
            using basic_t = typename std::remove_const<typename std::remove_reference<decltype(f)>::type>::type;
            *const_cast<basic_t *>(&f) = lexical_cast<basic_t>(vsl[i], status);
            if(!status)
              {
                ERROR("Failed to parse token: " << vsl[i] << " for parameter: " << this->m_name);
              }
          });
        }
      else
        {
          if(vsl.size() != 1)
            {
              ERROR("Constructing parameter: " << m_name << " only one argument expected, but " << vsl.size() - 1 << " found");
              for( auto const& v : vsl ) ERROR( v );
              return false;
            }
          m_value = lexical_cast<value_t>(vsl[0], status);
          if(status == false)
            {
              ERROR("Failed to parse token: " << vsl[0] << " for parameter: " << m_name);
              return false;
            }
        }
      return status;
    }

  private:
    void help(const value_t &def)
    {
      std::string type = type_string<value_t>();
      type = replaceAll(type, "AmpGen::", ""); /// remove namespaces
      if( std::is_same_v<value_t, std::string> ) 
        type = "string";
      if( std::is_same_v<value_t, std::vector<std::string>> ) 
        type = "strings";
      std::cout << " " << bold_on << std::left << std::setw(27) << m_name << bold_off << std::setw(20) << "[" + type + "]";
      auto tokens = split(m_helpString, '\n');
      if(tokens.size() == 0)
        std::cout << std::endl;
      for(size_t i = 0; i < tokens.size(); ++i)
        {
          if(i == 0)
            {
              std::cout << tokens[i];
              if constexpr(isVector<value_t>::value)
                {
                  if(def != value_t())
                    std::cout << " (default = {" << vectorToString(def, ",") << "})";
                }
              else if constexpr(isTuple<value_t>::value)
                {
                  if(def != value_t())
                    std::cout << " (default = [" << tupleToString(def, ",") << "])";
                }
              else
                {
                  if(def != value_t())
                    std::cout << " (default = " << def << ")";
                }
              std::cout << std::endl;
            }
          else
            std::cout << std::string(48, ' ') << tokens[i] << std::endl;
        }
    }
    bool setFromOptionsParser()
    {
      auto parser = OptionsParser::getMe();
      auto line = parser->find(m_name);
      if(line == parser->end())
        return false;
      std::vector<std::string> vsl = line->second;
      vsl.erase(vsl.begin());
      return setFromStrings(vsl);
    }
  };
  template <typename T> std::ostream &operator<<(std::ostream &os, const Property<T> &np);
  template <typename... T> std::string helpStringOptions(const std::string &header, const T &... args);
  template <typename T> bool operator==(const T &val, const Property<T> &prop) { return val == prop.value(); }
};

template <typename T> std::ostream &AmpGen::operator<<(std::ostream &os, const AmpGen::Property<T> &np)
{
  os << np.name() << " = ";
  if constexpr(!isVector<T>::value)
    os << np.m_value;
  else
    {
      for(unsigned i{0}; i != np.m_value.size(); ++i)
        os << np.m_value[i] << " ";
    }
  return os;
}

template <typename... T> std::string AmpGen::helpStringOptions(const std::string &header, const T &... args)
{
  std::stringstream rt;
  rt << header;
  for_each(std::make_tuple(args...), [&rt](const auto &f) mutable { rt << "\n\033[3m " << f.first << "\033[0m: " << f.second; });
  return rt.str();
}
#endif
