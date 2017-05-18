#ifndef ELEMENTCOUNTS_H_
#define ELEMENTCOUNTS_H_

#include "IDataStructure.h"
#include <vector>
#include <string>

class ElementCounts : public IDataStructure {
private:
    const std::vector<std::string> geneLabels;
    const std::function<std::string (unsigned)> elementLabelGenerator;
    unsigned curGene;

    std::unordered_map<unsigned, std::vector<int>> data;

public:
    ElementCounts(const std::function<std::string (unsigned)> elementLabelGenerator,
                  const std::vector<std::string> & geneLabels);

    virtual void sGeneInput(unsigned gene);
    virtual void sElementInput(unsigned element, int position);

    virtual const std::unordered_map<unsigned, std::vector<int> >& getStructure() const;

    virtual SEXP getSEXP() const;

    virtual unsigned getElementCount() const;
    virtual unsigned getGeneCount() const;
    virtual std::string getElementLabel(unsigned element) const;
    virtual const std::string& getGeneLabel(unsigned gene) const;
    virtual const std::vector<std::string> & getGeneLabels() const;
    virtual ~ElementCounts() {};
};

#endif /* ELEMENTCOUNTS_H_ */
