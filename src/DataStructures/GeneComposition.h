#ifndef GENECOMPOSITION_H_
#define GENECOMPOSITION_H_

#include "IDataStructure.h"
#include <vector>
#include <map>
#include <string>
#include <unordered_set>

class GeneComposition : public IDataStructure {
private:
    const std::vector<std::string> geneLabels;
    const std::function<std::string (unsigned)> elementLabelGenerator;
    std::unordered_set<unsigned> elements;
    int curGene, empty, k;

    std::unordered_map<unsigned, std::vector<int>> data;

public:
    GeneComposition(const std::function<std::string (unsigned)> elementLabelGenerator,
                   const std::vector<std::string> & geneLabels);

    virtual void sGeneInput(unsigned gene);
    virtual void sElementInput(unsigned element, int position);

    virtual const std::unordered_map<unsigned, std::vector<int>>& getStructure() const;
    virtual SEXP getSEXP() const;

    virtual unsigned getElementCount() const;
    virtual unsigned getGeneCount() const;
    virtual std::string getElementLabel(unsigned element) const;
    virtual const std::string& getGeneLabel(unsigned gene) const;
    virtual std::unordered_map<unsigned, std::string> getElementLabels() const;
    virtual const std::vector<std::string> & getGeneLabels() const;
    virtual ~GeneComposition() {};
};

#endif /* GENECOMPOSITION_H_ */
