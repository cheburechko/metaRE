#ifndef IDATASTRUCTURE_H_
#define IDATASTRUCTURE_H_

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

class IDataStructure {

public:
    virtual void sGeneInput(unsigned gene) = 0;
    virtual void sElementInput(unsigned element, int position) = 0;

    virtual const std::unordered_map<unsigned, std::vector<int> >& getStructure() const = 0;
    virtual SEXP getSEXP() const = 0;

    virtual unsigned getElementCount() const = 0;
    virtual unsigned getGeneCount() const = 0;
    virtual std::string getElementLabel(unsigned element) const = 0;
    virtual const std::string& getGeneLabel(unsigned gene) const = 0;

    virtual std::unordered_map<unsigned, std::string> getElementLabels() const {
        std::unordered_map<unsigned, std::string> result;
        for (auto it : getStructure()) {
            result[it.first] = getElementLabel(it.first);
        }
        return result;
    }

    virtual const std::vector<std::string> & getGeneLabels() const = 0;
    virtual ~IDataStructure() {};
};

#endif /* IDATASTRUCTURE_H_ */
