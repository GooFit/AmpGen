#ifndef AMPGEN_MINUITPARAMETERSET_H
#define AMPGEN_MINUITPARAMETERSET_H

#include <iostream>
#include <map>
#include <vector>

#include "AmpGen/MinuitParameter.h"

namespace AmpGen {
  class MinuitExpression;

  class MinuitParameterSet {
  public:
    typedef std::vector<MinuitParameter *>::iterator iterator;
    typedef std::vector<MinuitParameter *>::const_iterator const_iterator;

    MinuitParameterSet();
    explicit MinuitParameterSet(const std::vector<MinuitParameter *> &params);
    MinuitParameterSet(const MinuitParameterSet &other) = delete;
    ~MinuitParameterSet();
    MinuitParameterSet *clone() const;

    bool add(MinuitParameter *parPtr);
    MinuitParameter *add(const std::string &name, const Flag &flag, const double &mean, const double &sigma, const double &min = 0, const double &max = 0);
    bool unregister(MinuitParameter *patPtr);
    MinuitProxy addOrGet(const std::string &name, const Flag &flag, const double &mean, const double &sigma, const double &min = 0, const double &max = 0);
    void loadFromStream();
    void loadFromFile(const std::string &name);
    void set(const double *x, const std::vector<unsigned> &mapping, const double *errX = 0);
    void resetToInit();
    void print(std::ostream &os = std::cout) const;
    void printVariable(std::ostream &os = std::cout) const;
    void set(const MinuitParameterSet &mps);
    bool rename(const std::string &name, const std::string &new_name);
    unsigned int size() const;

    const_iterator cbegin() const;
    const_iterator cend() const;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    MinuitProxy at(const std::string &key);
    MinuitProxy at(const size_t &index) const;
    MinuitProxy operator[](const std::string &key);
    MinuitProxy operator[](const std::string &key) const;
    MinuitProxy operator[](const size_t &key);
    MinuitProxy find(const std::string &key) const;
    bool contains(const std::string &key) const;
    double operator()(const std::string &name);

    void setFromMinuit(const double *x);
    void setMapping(const std::vector<unsigned> &m);
    void setFromMinuitIndex(const unsigned index, double v);
    double getFromMinuitIndex(const unsigned index);

  private:
    void tryParameter(const std::vector<std::string> &line);
    void tryAlias(const std::vector<std::string> &line);
    bool addToEnd(MinuitParameter *parPtr);

    std::vector<MinuitParameter *> m_parameters;
    std::vector<unsigned> m_mapping;
    std::map<std::string, MinuitParameter *> m_keyAccess;
  };
} // namespace AmpGen
#endif
//
